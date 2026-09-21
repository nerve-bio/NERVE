#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('a')
parser.add_argument('b')
parser.add_argument('--tol', type=float, default=0.0)
args = parser.parse_args()

a_root, b_root = Path(args.a), Path(args.b)
if not a_root.is_dir() or not b_root.is_dir():
    raise SystemExit('Both epitope output directories must exist.')
a_files = {p.relative_to(a_root) for p in a_root.rglob('*.csv')}
b_files = {p.relative_to(b_root) for p in b_root.rglob('*.csv')}
if not a_files or not b_files:
    raise SystemExit('Both epitope output directories must contain CSV files.')
different = a_files != b_files

for relative in sorted(a_files & b_files):
    a = pd.read_csv(a_root / relative)
    b = pd.read_csv(b_root / relative)
    try:
        pd.testing.assert_frame_equal(a, b, check_exact=args.tol == 0,
                                      atol=args.tol, rtol=0)
    except AssertionError as error:
        different = True
        print(f'{relative}: DIFFERENT')
        print(str(error).splitlines()[0])

if a_files - b_files:
    print('Only in A:', *sorted(a_files - b_files), sep='\n  ')
if b_files - a_files:
    print('Only in B:', *sorted(b_files - a_files), sep='\n  ')
print('Epitope CSVs differ.' if different else
      f'Epitope CSVs identical within tol={args.tol:g}.')
raise SystemExit(1 if different else 0)
