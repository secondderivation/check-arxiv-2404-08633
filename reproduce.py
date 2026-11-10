"""Verify immutable inputs and re-run checks in a separate, retained work directory."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
from compare_results import compare_run

ROOT = Path(__file__).resolve().parent

def manifest_files(root):
    entries = []
    for line in (root / 'MANIFEST.txt').read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        digest, name = line.split('  ', 1)
        p = Path(name)
        if len(digest) != 64 or p.is_absolute() or '..' in p.parts:
            raise ValueError('Unsafe manifest entry')
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Integrity failure: {name}')
        entries.append(name)
    if not entries or len(entries) != len(set(entries)):
        raise ValueError('Empty or duplicate manifest')
    return entries

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--all', action='store_true', help='include all CPL sampler checks')
    args = ap.parse_args()
    files = manifest_files(ROOT)
    manifest_hash = hashlib.sha256((ROOT / 'MANIFEST.txt').read_bytes()).hexdigest()
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    out = ROOT / '_runs' / stamp
    work = out / 'work'
    for name in files + ['MANIFEST.txt']:
        dest = work / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, dest)
    if (ROOT / 'data').exists():
        shutil.copytree(ROOT / 'data', work / 'data')
    report = {'started_utc': stamp, 'complete_suite': args.all, 'manifest_sha256': manifest_hash,
              'python': sys.version, 'checks': []}
    started = time.monotonic()
    bad = False
    try:
        fetch = subprocess.run([sys.executable, str(work / 'fetch_data.py')], capture_output=True, text=True)
        (out / 'fetch.stdout').write_text(fetch.stdout)
        (out / 'fetch.stderr').write_text(fetch.stderr)
        if fetch.returncode:
            raise RuntimeError('Pinned data fetch failed; see fetch.stderr')
        runs = json.loads((work / 'runs.json').read_text())
        for i, run in enumerate(runs):
            if run.get('slow') and not args.all:
                report['checks'].append({'output': run['out'], 'status': 'SKIP'})
                print('SKIP', run['out'], '(use --all)', flush=True)
                continue
            cmd = shlex.split(run['cmd'])
            if cmd[0] != 'python3':
                raise ValueError('Only Python check commands are supported')
            cmd[0] = sys.executable
            t = time.monotonic()
            proc = subprocess.run(cmd, cwd=work / 'checks', text=True, capture_output=True)
            (out / f'{i:02d}.stdout').write_text(proc.stdout)
            (out / f'{i:02d}.stderr').write_text(proc.stderr)
            errors = [f'Exit {proc.returncode}'] if proc.returncode else compare_run(run, work, proc.stdout)
            status = 'FAIL' if errors else 'PASS'
            bad |= bool(errors)
            report['checks'].append({'output': run['out'], 'status': status, 'errors': errors,
                                      'seconds': time.monotonic() - t})
            print(status, run['out'], '; '.join(errors), flush=True)
    except Exception as exc:
        bad = True
        report['error'] = str(exc)
        print('FAIL', exc, flush=True)
    finally:
        try:
            manifest_files(ROOT)
            manifest_files(work)
            if hashlib.sha256((ROOT / 'MANIFEST.txt').read_bytes()).hexdigest() != manifest_hash:
                raise ValueError('Manifest changed during execution')
            report['integrity_after'] = 'PASS'
        except Exception as exc:
            bad = True
            report['integrity_after'] = str(exc)
        report['seconds'] = time.monotonic() - started
        report['status'] = 'FAIL' if bad else ('PASS' if args.all else 'PASS_CORE_ONLY')
        (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
        print('Evidence:', out, flush=True)
    return int(bad)

if __name__ == '__main__':
    sys.exit(main())
