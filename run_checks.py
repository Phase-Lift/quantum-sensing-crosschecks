#!/usr/bin/env python3
"""Run packaged benchmarks locally; never log in, commit, or upload anything."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true', help='Use original repetition counts.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    mode = 'full' if args.full else 'quick'
    output = root / 'outputs' / mode
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1')
    n = '20000' if args.full else '100'
    m = '5000' if args.full else '100'
    stages = root / 'benchmarks'
    jobs = [
        ('unit tests', stages / '03_dynamical_crosscheck',
         [sys.executable, '-m', 'unittest', '-v', 'test_crosscheck.py']),
        ('joint visibility', stages / '01_joint_visibility',
         [sys.executable, 'joint_visibility.py', '--trials', n,
          '--output', str(output / 'visibility_results.json')]),
        ('commutator controls', stages / '02_commutator_controls',
         [sys.executable, 'commutator_controls.py', '--repetitions', n,
          '--output', str(output / 'commutator')]),
        ('dynamical cross-check', stages / '03_dynamical_crosscheck',
         [sys.executable, 'oscillator_crosscheck.py', '--repetitions', m,
          '--output', str(output / 'crosscheck')]),
    ]
    for label, cwd, command in jobs:
        print(f'Running {label} ({mode})', flush=True)
        subprocess.run(command, cwd=cwd, env=env, check=True)
    print(f'All requested commands completed. Results are in outputs/{mode}/.')
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f'A benchmark command failed with exit code {exc.returncode}.', file=sys.stderr)
        raise SystemExit(exc.returncode or 1)
