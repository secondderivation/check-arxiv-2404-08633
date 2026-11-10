"""C7: Table II binned fits and the bin-to-bin shifts. Usage: python3 c7_binned.py [dr1|dr2]"""
import sys, json
import numpy as np
import bao, fit_grid as fg

release = sys.argv[1] if len(sys.argv) > 1 else 'dr1'
pv = json.load(open('../paper_values.json'))
paper = pv['table2_binned']
bins = [(0.1, 0.6), (0.6, 1.1), (1.1, 4.16)]
om_grid = np.linspace(0.0005, 0.9995, 2000)
h_grid = np.linspace(6000, 14000, 3201)


def shift(a, b):
    lo, hi = (a, b) if a[0] < b[0] else (b, a)
    return (hi[0] - lo[0]) / np.hypot(lo[1], hi[2])


for rz, tag in (((2, '2-dp z'), (None, '3-dp z')) if release == 'dr1' else ((None, '3-dp z'),)):
    rows, cov = bao.load(release, round_z=rz)
    out = []
    for k, (a, b) in enumerate(bins):
        rs, cv = bao.subset(rows, cov, lambda r: a < r[0] < b)
        s = fg.summarize(om_grid, h_grid, fg.posterior(rs, cv, om_grid, h_grid))
        out.append(s['om'])
        line = (f'C7 {release} [{tag}] {a}<z<{b} ({", ".join(f"{r[0]}:{r[2][:2]}" for r in rs)}): '
                f'Om = {s["om"][0]:.4f} +{s["om"][1]:.4f} -{s["om"][2]:.4f}, H0rd/100 = {s["h"][0]/100:.1f} +{s["h"][1]/100:.1f} -{s["h"][2]/100:.1f}')
        if release == 'dr1':
            c, ep, em = paper['rows'][k]['Om']; m = s['om'][0]
            tol = 0.25 * (ep if m > c else em)
            line += f'  | paper {c} +{ep} -{em}: Δ={m-c:+.4f} (tol {tol:.4f}) {"PASS" if abs(m-c) <= tol else "FAIL"}'
        print(line)
    s12, s23 = shift(out[0], out[1]), shift(out[1], out[2])
    msg = f'C7 {release} [{tag}] shifts: bin1->bin2 {s12:.2f} sigma, bin2->bin3 {s23:.2f} sigma (direction {"up" if out[2][0] > out[1][0] else "down"})'
    if release == 'dr1':
        p = paper['shifts_sigma']
        msg += (f' | paper {p["bin1_to_bin2"]}, {p["bin2_to_bin3"]}: '
                f'{"PASS" if abs(s12 - p["bin1_to_bin2"]) <= 0.2 and abs(s23 - p["bin2_to_bin3"]) <= 0.2 else "FAIL"}')
    print(msg)
