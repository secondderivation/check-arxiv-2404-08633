"""C9: within-DESI test of constant Ωm. Goodness of fit of one flat-ΛCDM (Ωm, H0rd) to all bins, and a
same-data nested tests against
Ωm varying between the Table II redshift groups (addendum A1; the first pass advertised a per-bin statistic it did not
implement and used leave-one-out, which is not a likelihood-ratio test).
Usage: python3 c9_constancy.py [dr1|dr2]"""
import sys
import numpy as np
from scipy.stats import chi2, norm
import bao, fit_grid as fg
from scipy.optimize import minimize


def nested(rows, cov, bounds, random_starts=0, seed=0):
    """Δχ² of (a) one Ωm per group with a common H0r_d and (b) one (Ωm, H0r_d) per group, against one common pair."""
    icov = np.linalg.inv(cov); d = np.array([r[1] for r in rows])
    grp = np.array([0 if r[0] < 0.6 else 1 if r[0] < 1.1 else 2 for r in rows])
    zmax = np.array([max(r[0] for r in rows if (0 if r[0] < 0.6 else 1 if r[0] < 1.1 else 2) == g) for g in range(3)])
    unit = lambda om: fg.unit_model(rows, om)

    def chi(oms, hs):
        if any(o <= max(bounds[0], -1 / ((1 + zm) ** 3 - 1)) or o >= bounds[1] for o, zm in zip(oms, zmax)):
            return 1e30
        m = np.zeros(len(rows))
        for g in range(3):
            m[grp == g] = unit(oms[g])[grp == g] / hs[g]
        r = d - m
        return r @ icov @ r

    rng = np.random.default_rng(seed)
    checks = {}

    def best(f, x0s, name=None):
        res = [minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-9, fatol=1e-10, maxiter=20000)) for x0 in x0s]
        b = min(res, key=lambda r: r.fun)
        if random_starts and name:
            # optimiser check (addendum A2): random starts must not find a lower minimum than the designed starts
            n = len(x0s[0]) // 2 if name == 'b' else None
            rs = []
            for _ in range(random_starts):
                x = np.array(x0s[0], float)
                if name == 'null': x = [rng.uniform(0.15, 0.6), rng.uniform(0.9, 1.1)]
                elif name == 'a': x = list(rng.uniform(0.15, 0.6, 3)) + [rng.uniform(0.9, 1.1)]
                else: x = [v for _ in range(3) for v in (rng.uniform(0.15, 0.6), rng.uniform(0.9, 1.1))]
                rs.append(minimize(f, x, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=40000)).fun)
            rs = np.array(rs)
            checks[name] = (b.fun - rs.min(), int(np.sum(rs - b.fun < 1e-6)), len(rs))
        return b
    null = best(lambda p: chi([p[0]] * 3, [p[1] * 1e4] * 3), [[0.3, 1.0], [0.25, 1.05]], 'null')
    a = best(lambda p: chi(p[:3], [p[3] * 1e4] * 3), [[0.4, 0.25, 0.32, 1.02], [0.3, 0.3, 0.3, 1.0], [0.6, 0.3, 0.35, 1.0]], 'a')
    b = best(lambda p: chi(p[0::2], p[1::2] * 1e4), [[0.5, 0.95, 0.25, 1.03, 0.33, 1.0], [0.3, 1.0, 0.3, 1.0, 0.3, 1.0],
                                                     [0.7, 0.9, 0.25, 1.05, 0.35, 0.98]], 'b')
    fmt = lambda v: '[' + ', '.join(f'{x:.4f}' for x in v) + ']'
    return {'(a) 3 Om, common H0rd': (null.fun - a.fun, 2, f'Om={fmt(a.x[:3])}, H0rd={a.x[3]*1e4:.1f}'),
            '(b) 3 (Om, H0rd) pairs': (null.fun - b.fun, 4, f'Om={fmt(b.x[0::2])}, H0rd={fmt(b.x[1::2]*1e4)}'),
            'null chi2min': (null.fun, 0, f'Om={null.x[0]:.4f}, H0rd={null.x[1]*1e4:.1f}'),
            **({'optimiser checks': checks} if random_starts else {})}

if __name__ == '__main__':
    release = sys.argv[1] if len(sys.argv) > 1 else 'dr1'
    om_grid = np.linspace(0.20, 0.42, 1101); h_grid = np.linspace(8500, 11500, 1501)
    for rz, tag in (((2, '2-dp z'), (None, '3-dp z')) if release == 'dr1' else ((None, '3-dp z'),)):
        rows, cov = bao.load(release, round_z=rz)
        s = fg.summarize(om_grid, h_grid, fg.posterior(rows, cov, om_grid, h_grid))
        dof = len(rows) - 2
        pte = chi2.sf(s['chi2min'], dof)
        print(f'C9 {release} [{tag}] global flat-ΛCDM fit: chi2min={s["chi2min"]:.2f}, dof={dof}, PTE={pte:.3f} '
              f'-> {norm.isf(pte/2):.2f} sigma (two-sided equivalent)')
        # Addendum A1: same-data nested tests replace leave-one-out (which compared chi2 on different data).
        # Groups as Table II: z < 0.6, 0.6-1.1, > 1.1. Wilks asymptotics; groups chosen after seeing the data.
        for bounds, btag in (((0.001, 0.999), 'Om in (0,1)'), ((-0.25, 3.0), 'Om in (-0.25,3), E^2>0')):
            out = nested(rows, cov, bounds)
            for name, (d2, k, pars) in out.items():
                p = chi2.sf(d2, k)
                if k == 0:
                    print(f'   nested [{btag}] null: χ²min = {d2:.3f}, {pars}'); continue
                print(f'   nested [{btag}] {name}: Δχ² = {d2:.3f} for {k} extra parameters, p = {p:.3f} '
                      f'({norm.isf(p/2):.2f} sigma two-sided equivalent); fitted {pars}')
