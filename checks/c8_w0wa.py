"""C8: DESI-alone flat w0waCDM, priors as printed (l.267): w0∈[-3,1], wa∈[-3,2], w0+wa<0; Ωm∈[0,1], H0rd∈[0,20000].
All four parameters are sampled with emcee, as in the paper.
Addendum A1 (2026-10-07): emcee's own random state is seeded (the first pass seeded only the starts), 4 independently
started ensembles per case, rank-normalised split R-hat on one fixed walker per ensemble; summaries in out-r2/ (A1) and out-r3/ with archived chains (A2).
Usage: python3 c8_w0wa.py [dr1|dr2] [rz]   (rz=2 rounds z_eff to 2 decimals)"""
import sys, json, os, hashlib, platform
import numpy as np
import bao, mcmc_diag as md

release = sys.argv[1] if len(sys.argv) > 1 else 'dr1'
rz = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != 'none' else None
NZ = 64  # fixed Gauss-Legendre nodes; checked against quad below
STEPS, BURN = 20000, 5000
THIN = 25  # archived chains: every 25th step of every walker (addendum A2); fixed walker kept in full


def make_like(rows, cov):
    icov = np.linalg.inv(cov)
    d = np.array([r[1] for r in rows])
    zs = np.array(sorted({r[0] for r in rows}))
    x, w = np.polynomial.legendre.leggauss(NZ)
    idx = np.array([np.searchsorted(zs, r[0]) for r in rows]); kind = np.array([r[2] for r in rows])

    def E(z, om, w0, wa):  # om, w0, wa: (n, 1)
        return np.sqrt(om * (1 + z) ** 3 + (1 - om) * (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z)))

    def dists(om, w0, wa):  # vectorised: arrays of length n -> (n, rows)
        om, w0, wa = (np.atleast_1d(v)[:, None] for v in (om, w0, wa))
        zz = (0.5 * zs[:, None] * (x + 1)).ravel()[None, :]
        I = (1.0 / E(zz, om, w0, wa)).reshape(len(om), len(zs), NZ) @ w * 0.5 * zs
        dm, dh = I[:, idx], 1.0 / E(zs[idx][None, :], om, w0, wa)
        out = np.where(kind == 'DM_over_rs', dm, np.where(kind == 'DH_over_rs', dh, np.cbrt(zs[idx] * dm ** 2 * dh)))
        return bao.C_KMS * out

    def logp(p):
        om, w0, wa, h = p.T
        ok = (0 < om) & (om < 1) & (-3 < w0) & (w0 < 1) & (-3 < wa) & (wa < 2) & (w0 + wa < 0) & (0 < h) & (h < 20000)
        out = np.full(len(p), -np.inf)
        if ok.any():
            r = d - dists(om[ok], w0[ok], wa[ok]) / h[ok, None]
            v = -0.5 * np.einsum('ij,jk,ik->i', r, icov, r)
            out[ok] = np.where(np.isfinite(v), v, -np.inf)
        return out
    return logp, dists


def run(rows, cov, label, seed):
    logp, dists = make_like(rows, cov)
    errs = [np.max(np.abs(dists(*p)[0] / bao.predict(rows, 1.0, bao.E_cpl, *p) - 1))
            for p in ((0.3, -0.8, -0.5), (0.3, -0.55, -1.8), (0.1, -2.8, 1.7), (0.99, 0.9, -2.9))]
    start = lambda rng, n: np.column_stack([rng.uniform(0.2, 0.45, n), rng.uniform(-1.3, -0.4, n),
                                            rng.uniform(-2.0, 0.2, n), rng.uniform(9500, 10800, n)])
    keep = {}
    ch, diag = md.run_ensembles(logp, start, 4, STEPS, BURN, seed=seed, keep=keep)
    attempts = [dict(steps=STEPS, max_rhat=float(max(diag['rhat'])), min_steps_over_tau=float(min(diag['steps_over_tau'])))]
    if not (np.all(np.isfinite(diag['rhat'])) and max(diag['rhat']) < 1.01 and min(diag['steps_over_tau']) >= 50):
        keep = {}
        ch, diag = md.run_ensembles(logp, start, 4, 40000, BURN, seed=seed, keep=keep)
        attempts.append(dict(steps=40000, max_rhat=float(max(diag['rhat'])), min_steps_over_tau=float(min(diag['steps_over_tau']))))
    tag = f'{release}_{"3dp" if rz is None else "2dp"}_{"all" if "with " in label else "no_lrg051"}'
    os.makedirs('chains', exist_ok=True)
    path = f'chains/c8_{tag}.npz'
    allc = np.array(keep['chains'])  # (ensembles, steps, walkers, ndim)
    np.savez_compressed(path, thinned=allc[:, ::THIN].astype(np.float32), fixed_walker=allc[:, :, 0, :].astype(np.float32),
                        thin=THIN, burn=BURN, seed=seed, names=np.array(['Om', 'w0', 'wa', 'H0rd']))
    sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    rhat_ok = max(diag['rhat']) < 1.01
    names = ['Om', 'w0', 'wa', 'H0rd']
    q = {n: np.percentile(ch[:, i], [15.9, 50, 84.1]) for i, n in enumerate(names)}
    parts = [f'{n} = {m:.3f} +{h-m:.3f} -{m-l:.3f}' for n, (l, m, h) in q.items()]
    print(f'C8 {release} z{"3dp" if rz is None else "2dp"} [{label}] n={len(rows)} max quad-err={max(errs):.1e} '
          f'{diag["ensembles"]}x32 walkers, steps={diag["steps"]} burn={BURN} | ' + '; '.join(parts) +
          f' | P(w0>-1)={np.mean(ch[:,1] > -1):.3f} | P(wa<-2.8)={np.mean(ch[:,2] < -2.8):.3f}')
    print('   ' + md.fmt_diag(diag, names) + f' | max Rhat < 1.01: {"PASS" if rhat_ok else "MISS"}')
    print(f'   chains: {path} sha256 {sha}')
    summ = dict(quad_error=float(max(errs)), attempts=attempts, label=label, release=release, z='3dp' if rz is None else '2dp', seed=seed,
                quantiles={n: list(map(float, v)) for n, v in q.items()}, chains=path, chains_sha256=sha, rhat_below_1p01=bool(rhat_ok),
                environment=dict(python=platform.python_version(), platform=platform.platform(), numpy=np.__version__,
                                 emcee=md.emcee.__version__, scipy=__import__('scipy').__version__),
                diag={k: (np.asarray(v).tolist() if hasattr(v, '__len__') else v) for k, v in diag.items()})
    return ch, summ


rows, cov = bao.load(release, round_z=rz)
ch_all, s1 = run(rows, cov, 'with LRG z=0.51', 81000)
r2, c2 = bao.subset(rows, cov, lambda r: abs(r[0] - 0.51) > 0.005)
ch_no, s2 = run(r2, c2, 'without LRG z=0.51', 82000)
dw = np.median(ch_no[:, 1]) - np.median(ch_all[:, 1])
lo, hi = np.percentile(ch_no[:, 1], [15.9, 84.1])
print(f'C8 {release}: removing the z=0.51 bin moves median w0 by {dw:+.3f}; -1 inside 68% without it: {lo <= -1 <= hi}')
os.makedirs('runtime-results', exist_ok=True)
json.dump([s1, s2], open(f'runtime-results/c8_{release}_{"3dp" if rz is None else "2dp"}.json', 'w'), indent=1)
# New chains are archived in chains/ (thinned; hashes in runtime-results/*.json). The same seeds reproduce them exactly only in the
# recorded environment: Codex's Mac run with the same seeds gave w0 medians up to 0.008 away.

if not all(s["rhat_below_1p01"] and min(s["diag"]["steps_over_tau"]) >= 50 for s in (s1, s2)):
    sys.exit("Convergence gate failed after the permitted extension")
