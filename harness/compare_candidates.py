#!/usr/bin/python3
"""
    Compare two vaccine_candidates.csv files as biology rather than as bytes.

    compare_states.py answers "did anything change"; this answers "does it
    matter". Two runs can differ in a score's last digits and still nominate the
    same proteins in the same order, which is a different situation from one that
    drops candidates.

    Usage:
        python3 harness/compare_candidates.py runs/a/vaccine_candidates.csv \\
                                              runs/b/vaccine_candidates.csv

    Exit status is 0 when the candidate sets match and every shared score agrees
    within --tol, 1 otherwise.
"""

import argparse
import os
import sys

import pandas as pd

ID_COLUMNS = ('uniprot_accession_code', 'id ', 'id')


def identifier_column(frame, path):
    """Pick the column identifying a candidate, preferring the accession."""
    for name in ID_COLUMNS:
        if name in frame.columns and frame[name].notna().any():
            return name
    raise SystemExit(f'{path}: no usable identifier column; found {list(frame.columns)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('a')
    parser.add_argument('b')
    parser.add_argument('--tol', type=float, default=0.0,
                        help='score tolerance (default 0 = exact)')
    args = parser.parse_args()

    a = pd.read_csv(args.a)
    b = pd.read_csv(args.b)
    label_a, label_b = os.path.basename(os.path.dirname(args.a)) or 'A', \
                       os.path.basename(os.path.dirname(args.b)) or 'B'

    column = identifier_column(a, args.a)
    a[column] = a[column].astype(str)
    b[column] = b[column].astype(str)
    set_a, set_b = set(a[column]), set(b[column])
    shared = set_a & set_b
    union = set_a | set_b

    print('=' * 74)
    print(f'CANDIDATE SETS   {label_a}: {len(set_a)}    {label_b}: {len(set_b)}')
    print('=' * 74)
    jaccard = (len(shared) / len(union)) if union else 1.0
    print(f'  shared={len(shared)}  jaccard={jaccard:.4f}')

    only_a, only_b = sorted(set_a - set_b), sorted(set_b - set_a)
    for label, missing in ((label_a, only_a), (label_b, only_b)):
        if missing:
            print(f'\n  ONLY in {label} ({len(missing)}):')
            for key in missing[:20]:
                print(f'    {key}')
            if len(missing) > 20:
                print(f'    ... and {len(missing) - 20} more')

    scores_agree = True
    if shared and 'score' in a.columns and 'score' in b.columns:
        ordered = sorted(shared)
        sa = a.drop_duplicates(column).set_index(column).loc[ordered, 'score']
        sb = b.drop_duplicates(column).set_index(column).loc[ordered, 'score']
        delta = (sa - sb).abs()
        over = int((delta > args.tol).sum())
        scores_agree = over == 0
        print('\n' + '=' * 74)
        print('SCORES AND RANKING (shared candidates)')
        print('=' * 74)
        print(f'  max|d|={delta.max():.3g}  mean|d|={delta.mean():.3g}  #(>{args.tol:g})={over}')
        if len(ordered) > 2:
            print(f'  correlation: pearson={sa.corr(sb):.6f}  spearman={sa.corr(sb, method="spearman"):.6f}')
        for n in (10, 25, 50, 100):
            if n <= min(len(a), len(b)):
                top_a = set(a.nlargest(n, 'score')[column])
                top_b = set(b.nlargest(n, 'score')[column])
                print(f'  top-{n:<4} overlap: {len(top_a & top_b)}/{n}')

    identical = (not only_a) and (not only_b) and scores_agree
    print()
    print('EQUIVALENT' if identical else 'DIFFERENT')
    return 0 if identical else 1


if __name__ == '__main__':
    sys.exit(main())
