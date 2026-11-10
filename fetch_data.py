"""Download the data files this package runs on from their pinned links and check their SHA-256."""
import hashlib, json, os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
for d in json.load(open(os.path.join(HERE, 'datasets.json')))['datasets']:
    for f in d.get('files', []):
        dest = os.path.join(HERE, f['local'])
        if not (os.path.exists(dest) and hashlib.sha256(open(dest, 'rb').read()).hexdigest() == f['sha256']):
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            urllib.request.urlretrieve(f['url'], dest)
        got = hashlib.sha256(open(dest, 'rb').read()).hexdigest()
        if got != f['sha256']:
            sys.exit(f"{f['local']}: SHA-256 {got} does not match {f['sha256']}")
        print(f"ok  {f['local']}")
