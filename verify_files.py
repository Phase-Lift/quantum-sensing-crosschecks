#!/usr/bin/env python3
"""Check release checksums. A match establishes integrity, not authenticity."""
from pathlib import Path
import hashlib
import sys

root = Path(__file__).resolve().parent
errors = []
manifest = root / 'SHA256SUMS.txt'
if not manifest.is_file():
    raise SystemExit('Missing SHA256SUMS.txt')
count = 0
for line in manifest.read_text(encoding='utf-8').splitlines():
    expected, name = line.split('  ', 1)
    file = (root / name).resolve()
    if not file.is_relative_to(root) or not file.is_file():
        errors.append(f'Missing or invalid path: {name}')
        continue
    actual = hashlib.sha256(file.read_bytes()).hexdigest()
    if actual != expected:
        errors.append(f'Checksum mismatch: {name}')
    count += 1
if errors:
    print('\n'.join(errors), file=sys.stderr)
    raise SystemExit(1)
print(f'PASS: {count} listed release files match their SHA-256 checksums.')
print('Additional unlisted files, identity metadata, and scientific validity are not checked.')
