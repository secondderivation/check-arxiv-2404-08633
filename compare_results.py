"""Release replay checks, distinct from the frozen scientific claim criteria.
See VALIDATION.md. Historical FAIL labels (e.g. C6) must reproduce unchanged.
"""
import json
import hashlib
import math
import re
from numbers import Real
import numpy as np

NUM = re.compile(r'(?<![\w.])[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')

def finite_numbers(value):
    """Reject NaN/infinity anywhere in structured evidence, including references."""
    if isinstance(value, dict):
        return all(finite_numbers(v) for v in value.values())
    if isinstance(value, (list, tuple, np.ndarray)):
        return all(finite_numbers(v) for v in value)
    return math.isfinite(value) if isinstance(value, Real) else True

def finite_vector(value, n):
    try:
        a = np.asarray(value, dtype=float)
        return a.shape == (n,) and bool(np.all(np.isfinite(a)))
    except (TypeError, ValueError):
        return False

def nonnegative_error(value):
    return isinstance(value, Real) and math.isfinite(value) and value >= 0

def text_compare(old, new, atol=1e-9, rtol=1e-9):
    if NUM.sub('#', old) != NUM.sub('#', new):
        return ['Output structure or verdict changed']
    a, b = [float(x) for x in NUM.findall(old)], [float(x) for x in NUM.findall(new)]
    if len(a) != len(b) or any(not math.isfinite(y) or not math.isclose(x, y, abs_tol=atol, rel_tol=rtol) for x, y in zip(a, b)):
        return ['Numeric replay tolerance exceeded']
    return []

def diag_errors(diag, n):
    errors = []
    for key, test in [('rhat', lambda x: 0 < x < 1.01), ('steps_over_tau', lambda x: x >= 50)]:
        vals = diag.get(key, [])
        if not finite_vector(vals, n) or any(not test(x) for x in vals):
            errors.append(f'Convergence gate: {key}')
    for key, test in [('tau', lambda x: x > 0), ('ensemble_median_spread', lambda x: x >= 0)]:
        if key in diag and (not finite_vector(diag[key], n) or any(not test(x) for x in diag[key])):
            errors.append(f'Invalid diagnostic: {key}')
    return errors

def c8_compare(old, new):
    errors = []
    if len(old) != 2 or len(new) != 2:
        return ['C8 requires both all-bin and removal cases']
    if not finite_numbers(old) or not finite_numbers(new):
        return ['C8 non-finite numeric evidence']
    for role, rows in [('reference', old), ('current', new)]:
        for row in rows:
            for name in ('Om', 'w0', 'wa', 'H0rd'):
                q = row.get('quantiles', {}).get(name)
                if not finite_vector(q, 3) or not np.all(np.diff(np.asarray(q, dtype=float)) > 0):
                    errors.append(f'C8 invalid {role} quantiles: {name}')
            if 'quad_error' in row and not nonnegative_error(row['quad_error']):
                errors.append(f'C8 invalid {role} quadrature error')
    if errors:
        return errors
    for ref, cur in zip(old, new):
        for key in ('label', 'release', 'z', 'seed'):
            if ref[key] != cur.get(key):
                errors.append(f'C8 identity mismatch: {key}')
        errors += diag_errors(cur['diag'], 4)
        if not nonnegative_error(cur.get('quad_error')) or cur['quad_error'] > 1e-8:
            errors.append('C8 quadrature control failed')
        for name in ('Om', 'w0', 'wa', 'H0rd'):
            a = np.asarray(ref['quantiles'][name], dtype=float)
            b = np.asarray(cur['quantiles'][name], dtype=float)
            if b.shape != (3,) or not np.all(np.isfinite(b)) or not np.all(np.diff(b) > 0):
                errors.append(f'C8 invalid quantiles: {name}')
                continue
            tol = .1 * min(a[1] - a[0], a[2] - a[1])
            if name == 'w0':
                tol = min(tol, .02)
            if np.max(np.abs(a - b)) > tol:
                errors.append(f'C8 {name} quantile replay exceeded {tol:.6g}')
    if new[0]['release'] == 'dr1' and new[0]['z'] == '2dp':
        a, b = [np.asarray(x['quantiles']['w0']) for x in new]
        if abs(a[1] + .56) > .10 or abs(b[1] + .98) > .15 or b[1] - a[1] > -.25 or not b[0] <= -1 <= b[2]:
            errors.append('C8 frozen paper criterion failed')
    return errors

def c4_compare(old, new, metrics):
    # The first 30 lines are deterministic prior tests. Subsequent sampler lines
    # are tested against exact grid posteriors, not an old random realization.
    errors = text_compare(old.split('C4-MCMC')[0], new.split('C4-MCMC')[0])
    old_mock = '\n'.join(x for x in old.splitlines() if x.startswith('   mock:'))
    new_mock = '\n'.join(x for x in new.splitlines() if x.startswith('   mock:'))
    errors += text_compare(old_mock, new_mock)
    if len(metrics) != 5 or [x['z'] for x in metrics] != [.51, .706, .93, 1.317, 2.33]:
        errors.append('C4 missing or mismatched bins')
    for m in metrics:
        if not finite_numbers(m):
            errors.append('C4 non-finite numeric evidence')
            continue
        errors += diag_errors(m['diag'], 2)
        if any(not finite_vector(m.get(k), 3) or not np.all(np.asarray(m[k], dtype=float)[1:] > 0)
               for k in ('grid', 'mcmc', 'paper')):
            errors.append('C4 invalid summary or interval offsets')
            continue
        grid, mc, paper = (np.asarray(m[k], dtype=float) for k in ('grid', 'mcmc', 'paper'))
        if not np.all(np.isfinite(mc)) or np.max(np.abs(grid - mc)) > .1 * min(grid[1:]):
            errors.append('C4 sampler/grid agreement failed')
        if abs(mc[0] - paper[0]) > .25 * (paper[1] if mc[0] > paper[0] else paper[2]) or np.any(np.abs(mc[1:] / paper[1:] - 1) > .2):
            errors.append('C4 frozen paper criterion failed')
        if (not nonnegative_error(m.get('quad_error')) or not nonnegative_error(m.get('grid_refinement'))
                or m['quad_error'] > 1e-8 or m['grid_refinement'] > .001):
            errors.append('C4 numerical control failed')
    # Preserve scientific PASS/FAIL words even on stochastic output.
    if re.findall(r'\b(?:PASS|FAIL|MISS)\b', old) != re.findall(r'\b(?:PASS|FAIL|MISS)\b', new):
        errors.append('C4 verdict changed')
    return errors

def c9_compare(old, new):
    errors = []
    oldlines, newlines = old.splitlines(), new.splitlines()
    if len(oldlines) != len(newlines):
        return ['C9 output shape changed']
    for a, b in zip(oldlines, newlines):
        # Only fitted nuisance displays have a units-aware tolerance; p/chi2 stay strict.
        marker = 'fitted Om='
        if marker in a and marker in b:
            ap, av = a.split(marker); bp, bv = b.split(marker)
            errors += text_compare(ap, bp)
            ao, ah = av.split(', H0rd='); bo, bh = bv.split(', H0rd=')
            errors += text_compare(ao, bo, atol=.0002, rtol=0)
            errors += text_compare(ah, bh, atol=.002, rtol=0)
        else:
            errors += text_compare(a, b)
    return errors

def compare_run(run, work, stdout):
    old = (work / 'checks' / run['out']).read_text()
    if 'c8_' in run['out']:
        name = run['out'].split('/')[-1].replace('.out', '.json')
        oldj = json.loads((work / 'checks' / 'out-r3' / name).read_text())
        newj = json.loads((work / 'checks' / 'runtime-results' / name).read_text())
        errors = c8_compare(oldj, newj)
        for row in newj:
            path = work / 'checks' / row['chains']
            if hashlib.sha256(path.read_bytes()).hexdigest() != row['chains_sha256']:
                errors.append('C8 runtime chain hash mismatch')
        return errors
    if 'c4.out' in run['out']:
        return c4_compare(old, stdout, json.loads((work / 'checks/runtime-results/c4.json').read_text()))
    if 'c9_' in run['out']:
        return c9_compare(old, stdout)
    return text_compare(old, stdout)
