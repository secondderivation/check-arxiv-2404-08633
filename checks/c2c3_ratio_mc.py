"""C2 (Table III) and C3 (p-value) by the paper's ratio method (l.225, l.245), on DESI DR1 or DR2 inputs.
Usage: python3 c2c3_ratio_mc.py [dr1|dr2]"""
import sys, json
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
from scipy.interpolate import interp1d
import bao

release = sys.argv[1] if len(sys.argv) > 1 else 'dr1'
rows, cov = bao.load(release)
pv = json.load(open('../paper_values.json'))
OM_HI = 50.0  # table end only; above it solve() brackets upward without a cap (addendum A1)


def solver(z):
    """Fast exact inverse of Eq. (2): dense monotonic table + brentq polish."""
    # start 1e-6 above the E(z)^2 = 0 edge: at 1e-9 the integrand is nearly singular and quad warned
    om_min = -1.0 / ((1 + z) ** 3 - 1) + 1e-6
    grid = np.concatenate([np.linspace(om_min, 3, 3000), np.linspace(3, OM_HI, 300)[1:]])
    vals = np.array([bao.ratio(z, o) for o in grid])
    assert np.all(np.diff(vals) > 0)
    inv = interp1d(vals, grid)

    def solve(t):
        if t <= vals[0]:
            return np.nan
        if t >= vals[-1]:
            # no finite cap: the ratio rises to a finite supremum as Om -> inf, so bracket upward
            hi = OM_HI
            while bao.ratio(z, hi) < t:
                hi *= 10
                if hi > 1e12:
                    return np.inf  # guard exhausted at 1e12: a numerical limit, not a proof that no root exists (A2)
            solve.above_table += 1
            return brentq(lambda o: bao.ratio(z, o) - t, OM_HI, hi, xtol=1e-13, rtol=1e-13)
        guess = float(inv(t))
        lo, hi = max(om_min, guess - 0.02), min(OM_HI, guess + 0.02)
        if (bao.ratio(z, lo) - t) * (bao.ratio(z, hi) - t) > 0:
            lo, hi = om_min, OM_HI
        return brentq(lambda o: bao.ratio(z, o) - t, lo, hi, xtol=1e-13)
    solve.above_table = 0
    solve.lo, solve.hi = om_min, np.inf
    return solve


def bin_draws(z, n, seed):
    i = [k for k, r in enumerate(rows) if abs(r[0] - z) < 1e-9 and r[2] == 'DM_over_rs'][0]
    j = [k for k, r in enumerate(rows) if abs(r[0] - z) < 1e-9 and r[2] == 'DH_over_rs'][0]
    mu = np.array([rows[i][1], rows[j][1]])
    C = cov[np.ix_([i, j], [i, j])]
    rng = np.random.default_rng(seed)
    x = rng.multivariate_normal(mu, C, size=n)
    return mu, C, x[:, 0] / x[:, 1]


aniso = sorted({r[0] for r in rows if r[2] == 'DM_over_rs'} & {r[0] for r in rows if r[2] == 'DH_over_rs'})
print(f'C2 release={release} anisotropic bins: {aniso}')
paper_rows = {round(r['z'], 2): r for r in pv['table3_ratio']['rows']}
results = {}
for k, z in enumerate(aniso):
    solve = solver(z)
    mu, C, t = bin_draws(z, 10000, seed=20261007 + k)
    om = np.array([solve(v) for v in t])
    bad = ~np.isfinite(om)
    good = om[~bad]
    resid = max(abs(bao.ratio(z, o) - v) for o, v in zip(good, t[~bad]))  # every draw (addendum A1)
    lo, med, hi = bao.pct(good)
    # promised alternative: clip no-root draws to the range edge instead of dropping them
    clipped = np.where(np.isnan(om), solve.lo, np.where(np.isinf(om), 1e12, om))
    clo, cmed, chi_ = bao.pct(clipped)
    # physical-domain sensitivities: drop or clip draws outside [0, 1]
    phys = good[(good >= 0) & (good <= 1)]
    plo, pmed, phi = bao.pct(phys)
    qlo, qmed, qhi = bao.pct(np.clip(good, 0, 1))
    r = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
    om_obs = solve(mu[0] / mu[1])
    line = (f'C2 z={z}: DM/rd={mu[0]:.3f}±{np.sqrt(C[0,0]):.3f} DH/rd={mu[1]:.3f}±{np.sqrt(C[1,1]):.3f} r={r:.3f} '
            f'| Om(at mean ratio)={om_obs:.4f} mean={good.mean():.4f} median={med:.4f} +{hi-med:.4f} -{med-lo:.4f} '
            f'| no-root draws {bad.sum()} (above Om=50: {solve.above_table}) | max reinsertion residual (all draws) {resid:.1e}'
            f'\n   clip-no-root median={cmed:.4f} +{chi_-cmed:.4f} -{cmed-clo:.4f}'
            f' | outside [0,1]: {np.mean((good<0)|(good>1)):.4f} | drop-outside median={pmed:.4f} +{phi-pmed:.4f} -{pmed-plo:.4f}'
            f' | clip-to-[0,1] median={qmed:.4f} +{qhi-qmed:.4f} -{qmed-qlo:.4f}')
    pr = paper_rows.get(round(z, 2))
    if pr and release == 'dr1':
        c, ep, em = pr['Om']
        dc = med - c
        tol = 0.25 * (ep if dc > 0 else em)
        ok_c = abs(dc) <= tol
        ok_e = abs((hi - med) / ep - 1) <= 0.2 and abs((med - lo) / em - 1) <= 0.2
        okm = abs(good.mean() - c) <= tol
        line += (f'\n   paper {c} +{ep} -{em}: median Δ={dc:+.4f} (tol ±{tol:.4f}) central {"PASS" if ok_c else "FAIL"}'
                 f' [mean Δ={good.mean()-c:+.4f} {"PASS" if okm else "FAIL"}]; errors ratio +{(hi-med)/ep:.2f} -{(med-lo)/em:.2f} {"PASS" if ok_e else "FAIL"}')
    print(line)
    results[z] = dict(median=med, lo=lo, hi=hi, mean=float(good.mean()))

# C3: 50,000 draws at the LRG bin z≈0.51
z = [a for a in aniso if abs(a - 0.51) < 0.01][0]
solve = solver(z)
for seed in (1, 2, 3):
    mu, C, t = bin_draws(z, 50000, seed=seed)
    om = np.array([solve(v) for v in t])
    p = np.mean(om < 0.322)
    se = np.sqrt(p * (1 - p) / len(om))
    print(f'C3 {release} seed={seed}: P(Om<0.322) = {p:.4f} ± {se:.4f} (binomial) '
          f'-> {norm.isf(p):.2f} sigma one-sided, {norm.isf(p/2):.2f} sigma two-sided; '
          f'P(Om<0.315) = {np.mean(om < 0.315):.4f}; no-root {(~np.isfinite(om)).sum()}')
# Analytic cross-check: the ratio itself is ~Gaussian; P(ratio < ratio(0.322)) by the delta method
i = [k for k, r in enumerate(rows) if abs(r[0] - z) < 1e-9 and r[2] == 'DM_over_rs'][0]
j = [k for k, r in enumerate(rows) if abs(r[0] - z) < 1e-9 and r[2] == 'DH_over_rs'][0]
m = np.array([rows[i][1], rows[j][1]]); C = cov[np.ix_([i, j], [i, j])]
g = np.array([1 / m[1], -m[0] / m[1] ** 2])
sr = np.sqrt(g @ C @ g)
thr = bao.ratio(z, 0.322)
pd = norm.cdf((thr - m[0] / m[1]) / sr)
print(f'C3 {release} delta-method: ratio={m[0]/m[1]:.4f}±{sr:.4f}, ratio(0.322)={thr:.4f}, P={pd:.4f} -> {norm.isf(pd):.2f} sigma one-sided')
# Exact for Gaussian (D_M, D_H) with D_H > 0: P(ratio < R) = P(D_M - R D_H < 0)   (addendum A1)
for th in (0.315, 0.322):
    R = bao.ratio(z, th); v = np.array([1.0, -R])
    zs = (v @ m) / np.sqrt(v @ C @ v); pe = norm.cdf(-zs)
    print(f'C3 {release} exact Gaussian tail P(Om<{th}) = {pe:.5f} -> {zs:.3f} sigma one-sided, {norm.isf(pe/2):.3f} two-sided'
          + (f'; rounds to {float(f"{zs:.1f}")} (criterion: 2.2) {"PASS" if abs(zs-2.2) < 0.05 + 1e-12 else "FAIL"}, p in [0.011, 0.019] {"PASS" if 0.011 <= pe <= 0.019 else "FAIL"}' if release == 'dr1' and th == 0.322 else ''))
