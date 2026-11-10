"""Negative controls for failures the original package silently missed."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from compare_results import text_compare, c9_compare, c8_compare, c4_compare, diag_errors
from reproduce import manifest_files

class ValidationTests(unittest.TestCase):
    def test_changed_verdict_rejected(self):
        self.assertTrue(text_compare('x=1 PASS', 'x=1 MISS'))
    def test_changed_value_rejected(self):
        self.assertTrue(text_compare('chi2=13.56', 'chi2=14.56'))
    def test_nonfinite_rejected(self):
        self.assertTrue(text_compare('x=1', 'x=1e999'))
        self.assertTrue(diag_errors({'rhat':[float('nan')], 'steps_over_tau':[100]}, 1))
    def test_fitted_display_tolerance_only(self):
        a='p = 0.537; fitted Om=[0.3161], H0rd=[9674.2726]'
        self.assertFalse(c9_compare(a,a.replace('9674.2726','9674.2727')))
        self.assertTrue(c9_compare(a,a.replace('0.537','0.538')))
        self.assertTrue(c9_compare(a,a.replace('9674.2726','9675.2726')))
    def pair(self):
        old=json.loads((Path(__file__).parent/'checks/out-r3/c8_dr1_2dp.json').read_text())
        new=copy.deepcopy(old)
        for row in new:
            row['quad_error']=1e-15
        return old,new
    def test_hashes_are_not_numeric_observables(self):
        old,new=self.pair()
        new[0]['chains_sha256']='0'*64
        self.assertFalse(c8_compare(old,new))
    def test_missing_case_rejected(self):
        old,new=self.pair()
        self.assertTrue(c8_compare(old,new[:1]))
    def test_posterior_and_convergence_changes_rejected(self):
        old,new=self.pair()
        new[0]['quantiles']['w0'][1]+=.1
        self.assertTrue(c8_compare(old,new))
        old,new=self.pair()
        new[0]['diag']['rhat'][0]=1.011
        self.assertTrue(c8_compare(old,new))
    def test_c8_nonfinite_nested_evidence_rejected(self):
        paths = [('quad_error',), ('seed',), ('diag', 'rhat', 0),
                 ('diag', 'steps_over_tau', 0), ('diag', 'tau', 0),
                 ('diag', 'ensemble_median_spread', 0)]
        paths += [('quantiles', name, i) for name in ('Om', 'w0', 'wa', 'H0rd') for i in range(3)]
        for side in ('reference', 'current'):
            for path in paths:
                for value in (float('nan'), float('inf'), float('-inf')):
                    with self.subTest(side=side, path=path, value=value):
                        old, new = self.pair()
                        row = old[0] if side == 'reference' else new[0]
                        target = row
                        for key in path[:-1]: target = target[key]
                        target[path[-1]] = value
                        self.assertTrue(c8_compare(old, new))
        old, new = self.pair()
        new[0]['attempts'] = [{'max_rhat': float('nan')}]
        self.assertTrue(c8_compare(old, new))
    def c4_metrics(self):
        return [dict(z=z, grid=[.3, .04, .03], mcmc=[.3, .04, .03], paper=[.3, .04, .03],
                     quad_error=1e-15, grid_refinement=.0001,
                     diag=dict(rhat=[1.001, 1.001], steps_over_tau=[100, 100],
                               tau=[20, 20], ensemble_median_spread=[.001, .001]),
                     attempts=[dict(steps=10000, max_rhat=1.02, min_steps_over_tau=100)])
                for z in [.51, .706, .93, 1.317, 2.33]]
    def test_c4_nonfinite_nested_evidence_rejected(self):
        paths = [('quad_error',), ('grid_refinement',), ('z',),
                 ('diag', 'rhat', 0), ('diag', 'steps_over_tau', 0),
                 ('diag', 'tau', 0), ('diag', 'ensemble_median_spread', 0),
                 ('attempts', 0, 'max_rhat')]
        paths += [(name, i) for name in ('grid', 'mcmc', 'paper') for i in range(3)]
        for path in paths:
            for value in (float('nan'), float('inf'), float('-inf')):
                with self.subTest(path=path, value=value):
                    metrics = self.c4_metrics()
                    target = metrics[0]
                    for key in path[:-1]: target = target[key]
                    target[path[-1]] = value
                    self.assertTrue(c4_compare('', '', metrics))
    def test_error_signs_and_interval_shapes(self):
        for key in ('quad_error', 'grid_refinement'):
            metrics = self.c4_metrics(); metrics[0][key] = -1e-15
            self.assertTrue(c4_compare('', '', metrics))
        for side in ('reference', 'current'):
            old, new = self.pair()
            (old if side == 'reference' else new)[0]['quad_error'] = -1e-15
            self.assertTrue(c8_compare(old, new))
            for q in ([0, 0, 1], [0, 1], [[0, 1, 2]]):
                old, new = self.pair()
                (old if side == 'reference' else new)[0]['quantiles']['w0'] = q
                self.assertTrue(c8_compare(old, new))
        for name in ('grid', 'mcmc', 'paper'):
            for q in ([.3, 0, .03], [.3, .04, -.03], [.3, .04], [[.3, .04, .03]]):
                metrics = self.c4_metrics(); metrics[0][name] = q
                self.assertTrue(c4_compare('', '', metrics))
        for name, value in [('tau', 0), ('ensemble_median_spread', -.001)]:
            metrics = self.c4_metrics(); metrics[0]['diag'][name][0] = value
            self.assertTrue(c4_compare('', '', metrics))
    def test_valid_structures_and_failed_attempt_history_preserved(self):
        old, new = self.pair()
        new[0]['attempts'] = [dict(steps=20000, max_rhat=1.02, min_steps_over_tau=40)]
        self.assertFalse(c8_compare(old, new))
        self.assertFalse(c4_compare('', '', self.c4_metrics()))
    def test_corrupted_manifest_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'sample').write_text('original')
            digest=hashlib.sha256(b'original').hexdigest()
            (p/'MANIFEST.txt').write_text(digest+'  sample\n')
            self.assertEqual(manifest_files(p),['sample'])
            (p/'sample').write_text('changed')
            with self.assertRaises(ValueError): manifest_files(p)
    def test_downloaded_data_hash_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'bad-data').write_bytes(b'wrong content')
            (root/'fetch_data.py').write_text((Path(__file__).parent/'fetch_data.py').read_text())
            dataset={'datasets':[{'files':[{'local':'data/sample', 'url':(root/'bad-data').as_uri(), 'sha256':'0'*64}]}]}
            (root/'datasets.json').write_text(json.dumps(dataset))
            proc=subprocess.run([sys.executable,str(root/'fetch_data.py')],capture_output=True,text=True)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn('does not match',proc.stderr)
    def test_manifest_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'MANIFEST.txt').write_text('0'*64+'  ../escape\n')
            with self.assertRaises(ValueError): manifest_files(p)

if __name__=='__main__': unittest.main()
