"""Shared model and data loading for the 2404.08633 checks. Written from the paper's Eqs. (1)-(2)
and conventions.json; no author code exists (G13)."""
import os
import numpy as np
from decimal import Decimal, ROUND_HALF_UP
from scipy.integrate import quad

C_KMS = 299792.458
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')


def load(release='dr1', round_z=None):
    """Return list of (z, value, kind) and the full covariance."""
    if release == 'dr1':
        base = os.path.join(DATA, 'desi-dr1', 'desi_2024_gaussian_bao_ALL_GCcomb_')
    else:
        base = os.path.join(DATA, 'desi-dr2', 'desi_gaussian_bao_ALL_GCcomb_')
    rows = []
    for line in open(base + 'mean.txt'):
        if line.startswith('#') or not line.strip():
            continue
        z, v, k = line.split()
        if round_z is not None:
            # half-up on the decimal string, as DESI v1 printed it (0.295 -> 0.30); float round() gives 0.29
            z = float(Decimal(z).quantize(Decimal(1).scaleb(-round_z), rounding=ROUND_HALF_UP))
        else:
            z = float(z)
        rows.append((z, float(v), k))
    cov = np.loadtxt(base + 'cov.txt')
    assert cov.shape == (len(rows), len(rows))
    return rows, cov


def E_lcdm(z, om):
    return np.sqrt(om * (1 + z) ** 3 + 1 - om)


def E_cpl(z, om, w0, wa):
    a = 1.0 / (1 + z)
    de = (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    return np.sqrt(om * (1 + z) ** 3 + (1 - om) * de)


def comoving_int(z, Efun, *p):
    """∫0^z dz'/E(z')"""
    return quad(lambda x: 1.0 / Efun(x, *p), 0, z, epsabs=1e-12, epsrel=1e-11)[0]


def ratio(z, om):
    """Eq. (2): (D_M/r_d)/(D_H/r_d) = E(z) ∫0^z dz'/E."""
    return E_lcdm(z, om) * comoving_int(z, E_lcdm, om)


def predict(rows, h0rd, Efun, *p):
    """Model D_X/r_d for each row; h0rd in km/s."""
    out = []
    cache = {}
    for z, _, k in rows:
        if z not in cache:
            I = comoving_int(z, Efun, *p)
            cache[z] = (C_KMS / h0rd * I, C_KMS / h0rd / Efun(z, *p))
        dm, dh = cache[z]
        if k == 'DM_over_rs':
            out.append(dm)
        elif k == 'DH_over_rs':
            out.append(dh)
        elif k == 'DV_over_rs':
            out.append((z * dm ** 2 * dh) ** (1 / 3))
        else:
            raise ValueError(k)
    return np.array(out)


def subset(rows, cov, keep):
    idx = [i for i, r in enumerate(rows) if keep(r)]
    return [rows[i] for i in idx], cov[np.ix_(idx, idx)]


def pct(samples, q=(15.9, 50, 84.1)):
    return np.percentile(samples, q)


def asym(lo, mid, hi):
    return mid, hi - mid, mid - lo
