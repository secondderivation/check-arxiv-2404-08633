"""C6: tension numbers, from the paper's printed values and from our C4/C5 grid fits (2-dp and 3-dp z_eff)."""
import json
import numpy as np
import bao, fit_grid as fg

pv = json.load(open('../paper_values.json'))
ext = pv['external_Om']; T = pv['tensions']


def tens(a, b):
    """a, b = (central, +err, -err); error on the facing side."""
    lo, hi = (a, b) if a[0] < b[0] else (b, a)
    return (hi[0] - lo[0]) / np.hypot(lo[1], hi[2])


def show(label, a, printed):
    for k, ref in printed.items():
        t = tens(a, ext[k])
        print(f'C6 {label} vs {k}: {t:.2f} sigma (paper ~{ref}) {"PASS" if abs(t - ref) <= (0.15 if ref == 2.0 else 0.1) else "FAIL"}')


lrg = pv['table1_mcmc']['rows'][0]['Om']
show('printed LRG0.51 Om', lrg, T['lrg051_vs'])
show('printed all-DESI Om', pv['full_sample']['all']['Om'], T['all_vs'])
show('printed no-LRG0.51 Om', pv['full_sample']['no_lrg051']['Om'], T['no_lrg051_vs'])

def verdicts(om, printed):
    out = []
    for k, ref in printed.items():
        t = tens(om, ext[k]); ok = abs(t - ref) <= (0.15 if ref == 2.0 else 0.1) + 1e-9
        out.append(f'{k} {t:.3f}σ vs {ref} (Δ {t-ref:+.3f}) {"PASS" if ok else "FAIL"}')
    return ', '.join(out)


# Refits judged by the same ±0.1 criterion (addendum A1). LRG row from our own C4 grid fit, prior Om in [0,1].
rows, cov = bao.load('dr1')
rs, cv = bao.subset(rows, cov, lambda r: abs(r[0] - 0.51) < 0.005 and r[2] != 'DV_over_rs')
g1, g2 = np.linspace(0.0005, 0.9995, 2000), np.linspace(3000, 20000, 3401)
lrg_fit = fg.summarize(g1, g2, fg.posterior(rs, cv, g1, g2))['om']
print(f'C6 refit LRG0.51 Om={lrg_fit[0]:.4f} +{lrg_fit[1]:.4f} -{lrg_fit[2]:.4f}: ' + verdicts(lrg_fit, T['lrg051_vs']))

om_grid = np.linspace(0.20, 0.42, 1101); h_grid = np.linspace(8500, 11500, 1501)
for rz, tag in ((2, '2-dp z'), (None, '3-dp z')):
    rows, cov = bao.load('dr1', round_z=rz)
    for drop in (False, True):
        rs, cv = (rows, cov) if not drop else bao.subset(rows, cov, lambda r: abs(r[0] - 0.51) > 0.005)
        om = fg.summarize(om_grid, h_grid, fg.posterior(rs, cv, om_grid, h_grid))['om']
        lab = f'ours {tag} {"no-LRG0.51" if drop else "all-DESI"} Om={om[0]:.4f}'
        print('C6 ' + lab + ': ' + ', '.join(f'{k} {tens(om, ext[k]):.2f}σ' for k in ext))
        if rz == 2:  # the printed numbers came from 2-dp z, so only these refits are judged
            print('   verdict: ' + verdicts(om, T['no_lrg051_vs' if drop else 'all_vs']))
