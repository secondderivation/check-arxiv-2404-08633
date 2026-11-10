"""C4: Table I per-bin MCMC numbers by grid posteriors; then the prior test.
Addendum A1 (2026-10-07): grid refinement (2x), a seeded multi-ensemble emcee cross-check with the printed priors,
and a noiseless mock recovery on the likelihood maximum."""
import json, os
import numpy as np
import bao, fit_grid as fg, mcmc_diag as md
from scipy.optimize import minimize

rows, cov = bao.load('dr1')
pv = json.load(open('../paper_values.json'))
paper = {round(r['z'], 2): r for r in pv['table1_mcmc']['rows']}
h_grid = np.linspace(3000, 20000, 3401)
aniso = sorted({r[0] for r in rows if r[2] == 'DM_over_rs'} & {r[0] for r in rows if r[2] == 'DH_over_rs'})


def check(name, val, ref):
    c, ep, em = ref; m, up, dn = val
    tol = 0.25 * (ep if m > c else em)
    okc = abs(m - c) <= tol
    oke = abs(up / ep - 1) <= 0.2 and abs(dn / em - 1) <= 0.2
    return (f'{name} ours {m:.4f} +{up:.4f} -{dn:.4f} vs paper {c} +{ep} -{em}: '
            f'Δ={m-c:+.4f} (tol {tol:.4f}) {"PASS" if okc else "FAIL"}; err ratios +{up/ep:.2f} -{dn/em:.2f} {"PASS" if oke else "FAIL"}')


for z in aniso:
    rs, cv = bao.subset(rows, cov, lambda r: abs(r[0] - z) < 1e-9)
    for label, om_grid in (('prior Om in [0,1]', np.linspace(0.0005, 0.9995, 1000)),
                           ('wide prior Om in (Om_min,3]', np.linspace(fg.om_min_for(rs) + 1e-4, 3.0, 3000))):
        s = fg.summarize(om_grid, h_grid, fg.posterior(rs, cv, om_grid, h_grid))
        h100 = tuple(x / 100 for x in s['h'])
        print(f'C4 z={z} [{label}] MLE Om={s["mle_om"]:.4f} chi2min={s["chi2min"]:.2e}')
        ref = paper.get(round(z, 2))
        print('   ' + check('Om', s['om'], ref['Om']))
        print('   ' + check('H0rd/100', h100, ref['H0rd100']))


# ---- addendum A1 ----
GLx, GLw = np.polynomial.legendre.leggauss(64)


def fast_model(z, om):
    """D_M/r_d, D_H/r_d at H0 r_d = 1 km/s, vectorised over om (Gauss-Legendre; checked against quad below)."""
    om = np.atleast_1d(om)[:, None]
    zz = 0.5 * z * (GLx + 1)
    E = lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om)
    dm = bao.C_KMS * 0.5 * z * np.sum(GLw / E(zz), axis=1)
    dh = bao.C_KMS / E(z)[:, 0]
    return dm, dh


metrics = []
for k, z in enumerate(aniso):
    rs, cv = bao.subset(rows, cov, lambda r: abs(r[0] - z) < 1e-9)
    kinds = [r[2] for r in rs]; d = np.array([r[1] for r in rs]); icov = np.linalg.inv(cv)
    qerr = max(abs(fast_model(z, o)[0][0] / bao.predict(rs, 1.0, bao.E_lcdm, o)[kinds.index('DM_over_rs')] - 1) for o in (0.05, 0.3, 0.95))
    # 1) grid refinement
    s1 = fg.summarize(np.linspace(0.0005, 0.9995, 1000), h_grid, fg.posterior(rs, cv, np.linspace(0.0005, 0.9995, 1000), h_grid))
    g2o, g2h = np.linspace(0.00025, 0.99975, 2000), np.linspace(3000, 20000, 6801)
    s2 = fg.summarize(g2o, g2h, fg.posterior(rs, cv, g2o, g2h))
    dref = np.max(np.abs(np.array(s1['om']) - np.array(s2['om'])))
    # 2) MCMC cross-check, priors as printed: Om in [0,1], H0rd in [0, 20000]

    def logp(p):
        om, h = p[:, 0], p[:, 1]
        ok = (om > 0) & (om < 1) & (h > 0) & (h < 20000)
        out = np.full(len(p), -np.inf)
        if ok.any():
            dm, dh = fast_model(z, om[ok])
            m = np.column_stack([dm if kk == 'DM_over_rs' else dh for kk in kinds]) / h[ok, None]
            r = d - m
            out[ok] = -0.5 * np.einsum('ij,jk,ik->i', r, icov, r)
        return out
    start = lambda rng, n: np.column_stack([rng.uniform(0.1, 0.9, n), rng.uniform(8000, 12000, n)])
    ch, diag = md.run_ensembles(logp, start, 2, steps=10000, burn=2000, seed=40000 + 10 * k)
    attempts = [dict(steps=10000, max_rhat=float(max(diag['rhat'])), min_steps_over_tau=float(min(diag['steps_over_tau'])))]
    if not (np.all(np.isfinite(diag['rhat'])) and max(diag['rhat']) < 1.01 and min(diag['steps_over_tau']) >= 50):
        ch, diag = md.run_ensembles(logp, start, 2, steps=20000, burn=2000, seed=40000 + 10 * k)
        attempts.append(dict(steps=20000, max_rhat=float(max(diag['rhat'])), min_steps_over_tau=float(min(diag['steps_over_tau']))))
    lo, me, hi = np.percentile(ch[:, 0], [15.9, 50, 84.1])
    ref = paper.get(round(z, 2))
    metrics.append(dict(attempts=attempts, z=z, grid=list(map(float, s2['om'])), mcmc=[float(me), float(hi-me), float(me-lo)], paper=ref['Om'], quad_error=float(qerr), grid_refinement=float(dref), diag={k: np.asarray(v).tolist() for k,v in diag.items()}))
    print(f'C4-MCMC z={z} quad-err={qerr:.1e} | grid Om {s2["om"][0]:.4f} +{s2["om"][1]:.4f} -{s2["om"][2]:.4f} '
          f'(1x vs 2x grid max diff {dref:.1e}) | MCMC Om {me:.4f} +{hi-me:.4f} -{me-lo:.4f} | MCMC-grid median diff {me-s2["om"][0]:+.4f}')
    print('   ' + check('MCMC Om', (me, hi - me, me - lo), ref['Om']))
    print('   ' + md.fmt_diag(diag, ['Om', 'H0rd']))
    # 3) noiseless mock at Om=0.3, H0rd=10000: likelihood maximum, optimiser started from the grid MLE
    mock_rows = [(r[0], v, r[2]) for r, v in zip(rs, bao.predict(rs, 10000.0, bao.E_lcdm, 0.3))]
    gm = fg.summarize(g2o, g2h, fg.posterior(mock_rows, cv, g2o, g2h))
    md_ = np.array([r[1] for r in mock_rows])
    chi = lambda p: (lambda r: r @ icov @ r)(md_ - bao.predict(mock_rows, p[1] * 1e4, bao.E_lcdm, p[0]))
    opt = minimize(chi, [gm['mle_om'], gm['mle_h'] / 1e4], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-20, maxiter=4000))
    print(f'   mock: grid MLE Om={gm["mle_om"]:.5f}; optimiser Om={opt.x[0]:.9f} H0rd={opt.x[1]*1e4:.4f} '
          f'|Om-0.3|={abs(opt.x[0]-0.3):.1e} {"PASS" if abs(opt.x[0]-0.3) < 1e-6 else "FAIL"} (tolerance 1e-6); posterior median {gm["om"][0]:.4f} (not a criterion)')

os.makedirs('runtime-results', exist_ok=True)
json.dump(metrics, open('runtime-results/c4.json', 'w'), indent=2)
