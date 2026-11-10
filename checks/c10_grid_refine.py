"""C10 (addendum A2, 2026-10-07): grid-refinement coverage that A1 promised but ran only for the baseline C4 grids.
Each grid used for a reported number is re-run at 2x resolution in both axes. Criterion: median and both errors move
by < 0.025 x the facing error (0.1 x the protocol tolerance). Also prints posterior mass near the grid edges.
The C9 nested optimisations are checked by the spread of chi2min over random starts (criterion < 1e-6)."""
import numpy as np
from scipy.optimize import minimize
import bao, fit_grid as fg


def double(g):
    return np.linspace(g[0], g[-1], 2 * len(g) - 1)


def judge(label, rows, cov, og, hg):
    a = fg.summarize(og, hg, fg.posterior(rows, cov, og, hg))
    b = fg.summarize(double(og), double(hg), fg.posterior(rows, cov, double(og), double(hg)))
    P = np.exp(fg.posterior(rows, cov, og, hg) - fg.posterior(rows, cov, og, hg).max())
    pom, ph = P.sum(1), P.sum(0)
    edge = max(pom[0], pom[-1]) / pom.max(), max(ph[0], ph[-1]) / ph.max()
    worst = 0.0
    for k in ('om', 'h'):
        m1, u1, d1 = a[k]; m2, u2, d2 = b[k]
        err = min(u2, d2)
        worst = max(worst, abs(m1 - m2) / err, abs(u1 - u2) / err, abs(d1 - d2) / err)
    print(f'C10 {label}: Om {a["om"][0]:.5f} +{a["om"][1]:.5f} -{a["om"][2]:.5f} -> 2x {b["om"][0]:.5f} +{b["om"][1]:.5f} -{b["om"][2]:.5f} '
          f'| worst shift {worst:.4f} x error (criterion < 0.025) {"PASS" if worst < 0.025 else "FAIL"} '
          f'| edge density Om {edge[0]:.1e}, H0rd {edge[1]:.1e}')


# C4 wide prior (the N1 run)
rows, cov = bao.load('dr1')
h4 = np.linspace(3000, 20000, 3401)
for z in (0.51, 0.706, 0.93, 1.317, 2.33):
    rs, cv = bao.subset(rows, cov, lambda r: abs(r[0] - z) < 1e-9)
    judge(f'C4 wide prior z={z}', rs, cv, np.linspace(fg.om_min_for(rs) + 1e-4, 3.0, 3000), h4)

# C5 / C6 / C9-global grids (same grid)
g5o, g5h = np.linspace(0.20, 0.42, 1101), np.linspace(8500, 11500, 1501)
for rel, rz, tag in (('dr1', 2, '2-dp'), ('dr1', None, '3-dp'), ('dr2', None, '3-dp')):
    rows, cov = bao.load(rel, round_z=rz)
    judge(f'C5 {rel} {tag} all', rows, cov, g5o, g5h)
    rs, cv = bao.subset(rows, cov, lambda r: abs(r[0] - 0.51) > 0.005)
    judge(f'C5 {rel} {tag} no z=0.51', rs, cv, g5o, g5h)

# C7 groups
g7o, g7h = np.linspace(0.0005, 0.9995, 2000), np.linspace(6000, 14000, 3201)
for rel, rz, tag in (('dr1', 2, '2-dp'), ('dr1', None, '3-dp'), ('dr2', None, '3-dp')):
    rows, cov = bao.load(rel, round_z=rz)
    for a_, b_ in ((0.1, 0.6), (0.6, 1.1), (1.1, 4.16)):
        rs, cv = bao.subset(rows, cov, lambda r: a_ < r[0] < b_)
        judge(f'C7 {rel} {tag} {a_}<z<{b_}', rs, cv, g7o, g7h)

# C9 nested optimisations: chi2min spread over 8 random starts per model
import c9_constancy as c9  # noqa: E402
for rel, rz, tag in (('dr1', 2, '2-dp'), ('dr1', None, '3-dp'), ('dr2', None, '3-dp')):
    rows, cov = bao.load(rel, round_z=rz)
    out = c9.nested(rows, cov, (0.001, 0.999), random_starts=8, seed=101)
    for name, (gain, n_hit, n) in out['optimiser checks'].items():
        ok = gain < 1e-6
        print(f'C10 C9 nested {rel} {tag} model {name}: best random-start chi2min is {gain:+.1e} below the reported one '
              f'({"PASS" if ok else "FAIL"}: criterion < 1e-6); {n_hit}/{n} random starts reach the reported minimum within 1e-6')
