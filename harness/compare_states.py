#!/usr/bin/python3
"""
    Diff two sets of per-stage snapshots, attribute by attribute.

    Snapshots are produced by running NERVE with --keep_snapshots True, which
    pickles the protein list after every stage into <working_dir>/snapshots/.
    Comparing them localises the first stage at which two runs diverge, which is
    far more informative than discovering that the final candidate lists differ.

    Usage:
        # whole run, stage by stage
        python3 harness/compare_states.py --dir runs/a/snapshots runs/b/snapshots
        # one stage, restricted to some attributes
        python3 harness/compare_states.py a/04_adhesin.pkl b/04_adhesin.pkl --attrs p_ad

    The default tolerance is 0, i.e. exact equality. That is deliberate: with the
    numeric thread pools pinned (NERVE.py does this at import) the pipeline is
    reproducible, so any difference at all is a real difference. A loose
    tolerance would hide exactly the small reorderings these comparisons exist to
    catch. Use --tol only to characterise a difference you have already found.

    Exit status is 0 when the two sides agree within tolerance, 1 otherwise.
"""

import argparse
import os
import pickle
import sys
import warnings

# Load repository classes used by snapshots.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_MISSING = object()


def load_proteins(path):
    """Load a snapshot pickle into a list of Protein objects."""
    with open(path, 'rb') as handle:
        # Saved NumPy pickles refer to a legacy module path. NumPy 2 warns
        # while resolving it; the loaded values are compared normally below.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning,
                                    message=r'numpy\.core\.numeric is deprecated.*')
            obj = pickle.load(handle)
    if isinstance(obj, dict) and 'list_of_proteins' in obj:
        return obj['list_of_proteins']
    return obj


def protein_key(protein):
    """Identify a protein across two runs."""
    return str(getattr(protein, 'accession', None) or getattr(protein, 'id', ''))


def keyed_proteins(path):
    result = {}
    for protein in load_proteins(path):
        key = protein_key(protein)
        if not key or key == 'None' or key in result:
            raise ValueError(f'{path}: missing or duplicate protein key: {key!r}')
        result[key] = attributes(protein)
    return result


def attributes(protein):
    """The comparable state of one protein."""
    return {k: v for k, v in vars(protein).items() if not k.startswith('_')}


def canonical_secondary_localization_ties(locations):
    """Ignore only the order of equally scored secondary PSORTb locations.

    The selected first location is used downstream and remains order-sensitive.
    Unequal-score ordering remains significant too. PSORTb has returned the
    same secondary scores in varying tie order across otherwise identical runs.
    """
    if not isinstance(locations, list) or len(locations) < 3:
        return locations
    if not all(hasattr(loc, 'localization') and hasattr(loc, 'reliability')
               for loc in locations):
        return locations
    tail = locations[1:]
    if any(left.reliability < right.reliability
           for left, right in zip(tail, tail[1:])):
        return locations
    return [locations[0], *sorted(tail, key=lambda loc: (-loc.reliability,
                                                          loc.localization))]


def values_equal(a, b, tol):
    """Compare two attribute values recursively.

    output: (equal, magnitude_of_difference_or_None, note)
    """
    if a is _MISSING or b is _MISSING:
        return False, None, 'present on one side only'
    if a is None or b is None:
        return (a is None and b is None), None, ('' if a is b else f'{a!r} vs {b!r}')
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b, None, ('' if a == b else f'{a!r} vs {b!r}')
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        magnitude = abs(float(a) - float(b))
        return magnitude <= tol, magnitude, ('' if magnitude <= tol else f'{a} vs {b} (d={magnitude:.3g})')
    if isinstance(a, str) and isinstance(b, str):
        return a == b, None, ('' if a == b else f'{a!r} vs {b!r}')
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return False, None, f'keys differ: {sorted(set(a) ^ set(b))}'
        worst = None
        for key in a:
            equal, magnitude, note = values_equal(a[key], b[key], tol)
            if not equal:
                return False, magnitude, f'[{key!r}] {note}'
            if magnitude is not None:
                worst = magnitude if worst is None else max(worst, magnitude)
        return True, worst, ''
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        # Preserve list order.
        if len(a) != len(b):
            return False, None, f'length {len(a)} vs {len(b)}'
        worst = None
        for position, (x, y) in enumerate(zip(a, b)):
            equal, magnitude, note = values_equal(x, y, tol)
            if not equal:
                return False, magnitude, f'[{position}] {note}'
            if magnitude is not None:
                worst = magnitude if worst is None else max(worst, magnitude)
        return True, worst, ''
    try:
        import numpy as np
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            first, second = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
            if first.shape != second.shape:
                return False, None, f'shape {first.shape} vs {second.shape}'
            if not first.size:
                return True, 0.0, ''
            magnitude = float(np.nanmax(np.abs(first - second)))
            return magnitude <= tol, magnitude, ('' if magnitude <= tol else f'max|d|={magnitude:.3g}')
    except ImportError:
        pass
    # Compare persisted object state.
    if hasattr(a, '__dict__') and hasattr(b, '__dict__'):
        if type(a).__module__ != type(b).__module__ or type(a).__name__ != type(b).__name__:
            return False, None, f'types differ: {type(a).__name__} vs {type(b).__name__}'
        return values_equal(attributes(a), attributes(b), tol)
    equal = (a == b)
    return equal, None, ('' if equal else f'{a!r} vs {b!r}')


def compare_snapshots(path_a, path_b, tol=0.0, attrs=None,
                      canonical_localization_ties=False):
    """Compare two snapshot files; returns a result dict."""
    a = keyed_proteins(path_a)
    b = keyed_proteins(path_b)
    common = sorted(set(a) & set(b))
    per_attr, examples = {}, {}
    for key in common:
        names = attrs if attrs else sorted(set(a[key]) | set(b[key]))
        for name in names:
            left, right = a[key].get(name, _MISSING), b[key].get(name, _MISSING)
            if canonical_localization_ties and name == 'localization':
                left = canonical_secondary_localization_ties(left)
                right = canonical_secondary_localization_ties(right)
            equal, magnitude, note = values_equal(left, right, tol)
            if equal:
                continue
            entry = per_attr.setdefault(name, {'count': 0, 'worst': None})
            entry['count'] += 1
            if magnitude is not None:
                entry['worst'] = magnitude if entry['worst'] is None else max(entry['worst'], magnitude)
            examples.setdefault(name, []).append((key, note))
    return {'only_a': sorted(set(a) - set(b)), 'only_b': sorted(set(b) - set(a)),
            'n_common': len(common), 'n_a': len(a), 'n_b': len(b),
            'per_attr': per_attr, 'examples': examples}


def print_result(result, label='', tol=0.0):
    """Print one comparison; returns True when the two sides agree."""
    head = f'{label}: ' if label else ''
    agrees = not result['per_attr'] and not result['only_a'] and not result['only_b']
    if agrees:
        print(f'{head}IDENTICAL  ({result["n_common"]} proteins, tol={tol:g})')
        return True
    print(f'{head}DIFFERENT  (A={result["n_a"]} B={result["n_b"]} '
          f'common={result["n_common"]}, tol={tol:g})')
    for side in ('only_a', 'only_b'):
        if result[side]:
            shown = ', '.join(result[side][:8])
            print(f'    {side}: {len(result[side])} ({shown}{"..." if len(result[side]) > 8 else ""})')
    for name in sorted(result['per_attr'], key=lambda n: -result['per_attr'][n]['count']):
        entry = result['per_attr'][name]
        worst = f'  max|d|={entry["worst"]:.3g}' if entry['worst'] is not None else ''
        print(f'    {name:<48} {entry["count"]:>6} proteins{worst}')
        for key, note in result['examples'][name][:2]:
            if note:
                print(f'        e.g. {key}: {note}')
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('a', help='snapshot .pkl, or directory with --dir')
    parser.add_argument('b', help='snapshot .pkl, or directory with --dir')
    parser.add_argument('--dir', action='store_true', help='compare whole snapshot directories')
    parser.add_argument('--tol', type=float, default=0.0,
                        help='numeric tolerance (default 0 = exact)')
    parser.add_argument('--attrs', default=None,
                        help='comma-separated attributes to restrict the comparison to')
    parser.add_argument('--canonical-localization-ties', action='store_true',
                        help='ignore equal-score order among secondary PSORTb locations')
    args = parser.parse_args()
    attrs = [a.strip() for a in args.attrs.split(',')] if args.attrs else None

    if not args.dir:
        result = compare_snapshots(args.a, args.b, args.tol, attrs,
                                   args.canonical_localization_ties)
        return 0 if print_result(result, os.path.basename(args.a), args.tol) else 1

    in_a = {f for f in os.listdir(args.a) if f.endswith('.pkl')}
    in_b = {f for f in os.listdir(args.b) if f.endswith('.pkl')}
    for name in sorted(in_a ^ in_b):
        print(f'{name}: present in only one run')
    agree = not (in_a ^ in_b)
    for name in sorted(in_a & in_b):
        result = compare_snapshots(os.path.join(args.a, name), os.path.join(args.b, name),
                                   args.tol, attrs, args.canonical_localization_ties)
        agree &= print_result(result, name, args.tol)
    print()
    print(f'{"All stages identical" if agree else "Stages differ"} (tol={args.tol:g})')
    return 0 if agree else 1


if __name__ == '__main__':
    sys.exit(main())
