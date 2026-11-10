"""Explicit evidence for the two minor prior/notation notes; no new adverse claim."""
import numpy as np
import bao
import fit_grid as fg

rows, cov = bao.load('dr1')
rows, cov = bao.subset(rows, cov, lambda r: abs(r[0] - .51) < .005)
h = np.linspace(3000, 20000, 3401)
for label, om in [('printed', np.linspace(.0005, .9995, 1000)),
                  ('wide', np.linspace(fg.om_min_for(rows) + .0001, 3, 3000))]:
    s = fg.summarize(om, h, fg.posterior(rows, cov, om, h))
    m, up, dn = s['om']
    sigma = (m - .315) / np.hypot(dn, .007)
    print(f'N1 {label} prior: Om={m:.4f} +{up:.4f} -{dn:.4f}; facing-error Planck discrepancy={sigma:.3f} sigma')
print('N4 definitions: alpha_parallel=(DH/rd)/(DH_fid/rd_fid); alpha_perp=(DM/rd)/(DM_fid/rd_fid)')
print('N4 DESI convention: alpha_AP=alpha_parallel/alpha_perp')
print('N4 therefore alpha_AP=(DH/DM)/(DH_fid/DM_fid), and DM/DH=(DM_fid/DH_fid)/alpha_AP')
print('N4 this is a notation clarification; no distance equation or numerical fit is changed.')
