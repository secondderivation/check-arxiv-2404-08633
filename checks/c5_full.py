"""C5: flat ΛCDM fit to all DESI DR1 BAO (full covariance), with 2- and 3-decimal z_eff, with and without z=0.51.
Usage: python3 c5_full.py [dr1|dr2]"""
import sys, json
import numpy as np
import bao, fit_grid as fg

release = sys.argv[1] if len(sys.argv) > 1 else 'dr1'
pv = json.load(open('../paper_values.json'))
om_grid = np.linspace(0.20, 0.42, 1101)
h_grid = np.linspace(8500, 11500, 1501)
res = {}
for rz in ((None, 2) if release == 'dr1' else (None,)):
    rows, cov = bao.load(release, round_z=rz)
    tag = '3-dp z' if rz is None else '2-dp z'
    for drop in (False, True):
        rs, cv = (rows, cov) if not drop else bao.subset(rows, cov, lambda r: abs(r[0] - 0.51) > 0.005)
        s = fg.summarize(om_grid, h_grid, fg.posterior(rs, cv, om_grid, h_grid))
        m, up, dn = s['om']
        res[(tag, drop)] = s
        print(f'C5 {release} [{tag}, {"without" if drop else "with"} z=0.51, n={len(rs)}] Om = {m:.4f} +{up:.4f} -{dn:.4f}, '
              f'MLE {s["mle_om"]:.4f}, H0rd = {s["h"][0]:.1f} +{s["h"][1]:.1f} -{s["h"][2]:.1f}, chi2min = {s["chi2min"]:.2f}')

if release == 'dr1':
    fs = pv['full_sample']
    for drop, key in ((False, 'all'), (True, 'no_lrg051')):
        c, ep, em = fs[key]['Om']; mle = fs[key]['mle']
        for tag in ('2-dp z', '3-dp z'):
            s = res[(tag, drop)]; m, up, dn = s['om']
            ok = abs(m - c) <= 0.003 and abs(up - ep) <= 0.002 and abs(dn - em) <= 0.002 and abs(s['mle_om'] - mle) <= 0.002
            print(f'C5 check {key} [{tag}]: ours {m:.4f} +{up:.4f} -{dn:.4f} MLE {s["mle_om"]:.4f} vs paper {c} +{ep} -{em} MLE {mle}: {"PASS" if ok else "FAIL"}')
    shift = res[('3-dp z', False)]['om'][0] - res[('2-dp z', False)]['om'][0]
    print(f'C5 footnote: Om(3-dp z) - Om(2-dp z) = {shift:+.4f}; paper needs +0.0025..+0.01 to explain 0.290 vs DESI 0.295: '
          f'{"PASS" if 0.0025 <= shift <= 0.01 else "FAIL"}')
