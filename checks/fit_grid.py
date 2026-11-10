"""Exact 2D posteriors (Ωm, H0 r_d) for flat ΛCDM on a dense grid; marginals by summation, quantiles from a trapezoid-integrated CDF.
Distances scale as 1/(H0 r_d), so for each Ωm the H0 r_d dependence is evaluated on a grid in u = 1/(H0 r_d)."""
import numpy as np
import bao


def unit_model(rows, om):
    """D_X/r_d at H0 r_d = 1 km/s (i.e. c * dimensionless distance)."""
    return bao.predict(rows, 1.0, bao.E_lcdm, om)


def posterior(rows, cov, om_grid, h_grid):
    icov = np.linalg.inv(cov)
    d = np.array([r[1] for r in rows])
    logL = np.empty((len(om_grid), len(h_grid)))
    for a, om in enumerate(om_grid):
        m1 = unit_model(rows, om)
        # model = m1 / h; chi2 = (d - m1/h)^T C^-1 (d - m1/h)
        A = m1 @ icov @ m1; B = m1 @ icov @ d; Cc = d @ icov @ d
        u = 1.0 / h_grid
        logL[a] = -0.5 * (Cc - 2 * B * u + A * u ** 2)
    return logL


def summarize(om_grid, h_grid, logL):
    P = np.exp(logL - logL.max())
    pom = P.sum(axis=1); ph = P.sum(axis=0)
    def ci(x, p):
        # Trapezoid-integrated CDF on the grid points (fixed 2026-10-07 after Codex B: a plain cumsum
        # interpolated at grid points displaced every quantile by about half a grid step).
        c = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(x))]); c /= c[-1]
        lo, med, hi = np.interp([0.159, 0.5, 0.841], c, x)
        return med, hi - med, med - lo
    a, b = np.unravel_index(np.argmax(logL), logL.shape)
    return dict(om=ci(om_grid, pom), h=ci(h_grid, ph), mle_om=om_grid[a], mle_h=h_grid[b], chi2min=-2 * logL.max())


def om_min_for(rows):
    zmax = max(r[0] for r in rows)
    return -1.0 / ((1 + zmax) ** 3 - 1)
