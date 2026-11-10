"""C1: re-derive Eq. (2) from Eq. (1), take limits, check monotonicity and the reachable ratio range."""
import sympy as sp
import numpy as np
from scipy.optimize import brentq
import bao

z, zp, om, H0, c, rd = sp.symbols('z zp Omega_m H0 c r_d', positive=True)
I = sp.Function('I')(z)  # stands for the unevaluated integral ∫0^z dz'/E(z')
E = lambda x: sp.sqrt(om * (1 + x) ** 3 + 1 - om)
DM = c / H0 * I
DH = c / (H0 * E(z))
r = sp.cancel((DM / rd) / (DH / rd))
print('C1 derived ratio (I = ∫0^z dz\'/E):', r)
print('C1 equals Eq.(2) printed form E(z)*I:', sp.simplify(r - E(z) * I) == 0)
print('C1 depends on H0, r_d, c:', any(s in r.free_symbols for s in (H0, rd, c)))

# limits (G9): substitute into the integrand first, then integrate
def ratio_at(omv):
    Ei = sp.sqrt(omv * (1 + zp) ** 3 + 1 - omv)
    return sp.simplify(sp.sqrt(omv * (1 + z) ** 3 + 1 - omv) * sp.integrate(1 / Ei, (zp, 0, z)))
r0 = ratio_at(0)
print('C1 limit Om=0:', r0, ' equals z:', sp.simplify(r0 - z) == 0)
r1 = ratio_at(1)
print('C1 limit Om=1:', r1, ' equals 2[(1+z)^{3/2}-(1+z)]:',
      sp.simplify(r1 - 2 * ((1 + z) ** sp.Rational(3, 2) - (1 + z))) == 0)
# small z: ratio = z + a z^2 + ...; derive a from the integrand expansion
Eser = sp.series(E(zp), zp, 0, 3).removeO()
Iser = sp.integrate(sp.series(1 / E(zp), zp, 0, 3).removeO(), (zp, 0, z))
rser = sp.expand(sp.series(Eser.subs(zp, z) * Iser, z, 0, 3).removeO())
print('C1 small-z series:', sp.collect(sp.simplify(rser), z), '-> ratio -> z as z -> 0, Om enters at O(z^2)')

# monotonicity: d ratio / d Om on a grid (numeric). E^2 > 0 on [0, z] needs Om > -1/((1+z)^3 - 1).
for zz in (0.51, 0.706, 0.93, 1.317, 2.33):
    om_min = -1.0 / ((1 + zz) ** 3 - 1)
    oms = np.linspace(om_min + 1e-6, 3.0, 600)
    vals = np.array([bao.ratio(zz, o) for o in oms])
    d = np.diff(vals)
    sign = 'increasing' if np.all(d > 0) else ('decreasing' if np.all(d < 0) else 'NOT monotonic')
    print(f'C1 z={zz}: Om > {om_min:.4f} keeps E^2>0; ratio is {sign} in Om on [{oms[0]:.4f}, 3]; '
          f'ratio(Om_min)={vals[0]:.5f}, ratio(0)={bao.ratio(zz, 0.0):.5f}, ratio(0.315)={bao.ratio(zz, 0.315):.5f}, '
          f'ratio(1)={bao.ratio(zz, 1.0):.5f}, ratio(3)={vals[-1]:.5f}; '
          f'd ratio/d Om at 0.315 = {(bao.ratio(zz, 0.316) - bao.ratio(zz, 0.314)) / 0.002:.4f}')

# known-answer control: invert noiseless ratio at Om=0.3
for zz in (0.51, 0.706, 0.93, 1.317, 2.33):
    t = bao.ratio(zz, 0.3)
    s = brentq(lambda o: bao.ratio(zz, o) - t, 0.0, 1.0, xtol=1e-14)
    print(f'C1 control z={zz}: solved Om={s:.10f} residual={bao.ratio(zz, s) - t:.2e}')
