import random
import sys
import types

import pandas as pd
import pytest

try:
    from code import Adhesin
except ModuleNotFoundError as error:
    if error.name != 'tensorflow':
        raise
    sys.modules['tensorflow'] = types.SimpleNamespace(keras=None)
    from code import Adhesin
from code.Epitope import _coverage
from code.Function import DeepFriParser
from code.Protein import Protein


def old_windows(match, query, size, substitutions, mismatches):
    if substitutions >= size or mismatches >= size:
        return []
    return [dict(match=match[i:i + size], query=query, start_pos=i)
            for i in range(len(match) - size + 1)
            if match[i:i + size].count('+') <= substitutions
            and match[i:i + size].count(' ') <= mismatches]


def old_coverage(binders, sequence, length):
    result = {}
    for allele in binders.allele.unique():
        coverage = [0 for _ in sequence]
        for start in binders[binders.allele == allele].pos.to_list():
            for position in range(start, min(start + length, len(sequence) + 1)):
                coverage[position - 1] += 1
        result[allele] = ' '.join(map(str, coverage))
    return result


def old_deepfri(path, threshold=0.5):
    frame = pd.read_csv(path, comment='#')
    rows = []
    for protein in frame['Protein'].unique():
        group = frame[frame['Protein'] == protein]
        entries = [(row['Score'], row['GO_term/EC_number name'])
                   for _, row in group.iterrows()]
        rows.append([protein, ' | '.join(go for score, go in
                     sorted(entries, key=lambda entry: entry[0], reverse=True)
                     if score > threshold)])
    return pd.DataFrame(rows, columns=['Protein', 'Function'])


def test_rolling_windows_match_old_parser():
    rng = random.Random(7)
    for size in (1, 2, 9, 15):
        for length in (0, size - 1, size, size + 1, 100):
            match = ''.join(rng.choice('|+ ABC') for _ in range(length))
            for substitutions, mismatches in ((0, 0), (1, 1), (3, 1)):
                if substitutions >= size or mismatches >= size:
                    continue
                assert Protein.hsp_match_parser(match, 'QUERY', size,
                                                substitutions, mismatches) == old_windows(
                    match, 'QUERY', size, substitutions, mismatches)
    assert Protein.hsp_match_parser('||||', 'QUERY', 5, 1, 1) == []


def test_coverage_matches_old_output():
    rng = random.Random(11)
    for size in (0, 1, 9, 100):
        for length in (1, 9, 15):
            positions = [rng.randrange(1, size + 4) for _ in range(100)]
            binders = pd.DataFrame({'allele': ['B', 'A'] * 50, 'pos': positions})
            assert _coverage(binders, 'A' * size, length) == old_coverage(
                binders, 'A' * size, length)
    binders = pd.DataFrame({'allele': ['A', 'A', 'B'], 'pos': [1, 1, 0]})
    assert _coverage(binders, 'AAAA', 2) == old_coverage(binders, 'AAAA', 2)
    with pytest.raises(IndexError):
        _coverage(pd.DataFrame({'allele': ['A'], 'pos': [-10]}), 'AAAA', 2)


def test_deepfri_parser_matches_old_order_and_ties(tmp_path):
    path = tmp_path / 'predictions.csv'
    frame = pd.DataFrame([
        ('B', 0.7, 'first'), ('A', 0.5, 'excluded'),
        ('B', 0.9, 'high'), ('A', 0.8, 'kept'), ('B', 0.7, 'second'),
    ], columns=['Protein', 'Score', 'GO_term/EC_number name'])
    frame.to_csv(path, index=False)
    pd.testing.assert_frame_equal(DeepFriParser(path), old_deepfri(path))
    assert list(DeepFriParser(path)['Protein']) == ['B', 'A']
    frame.loc[0, 'Protein'] = None
    frame.to_csv(path, index=False)
    with pytest.raises(IndexError):
        DeepFriParser(path)


def test_selected_id_set_preserves_discard_order():
    for selected in ([], ['B'], ['B', 'B']):
        all_ids = ['A', 'B', 'B', 'C']
        assert [item for item in all_ids if item not in set(selected)] == [
            item for item in all_ids if item not in selected]


def test_streamed_ifeature_rows_and_duplicate_rejection(tmp_path, monkeypatch):
    proteins = [Protein('sp|P1|P1', 'AAAA'), Protein('sp|P2|P2', 'CCCC')]
    features = ('AAC', 'DPC', 'CTDC', 'CTDT', 'CTDD')

    def write_outputs(_args):
        feature = _args[_args.index('--type') + 1]
        path = _args[_args.index('--out') + 1]
        with open(path, 'w') as output:
            output.write('header\n')
            output.write(f'P2\t{features.index(feature) + 2}.0\n')
            output.write(f'P1\t{features.index(feature) + 1}.0\n')

    monkeypatch.setattr(Adhesin, 'run_command', write_outputs)
    Adhesin.extract_features(proteins, '', '', str(tmp_path), '')
    assert list(proteins[0].model_raw_data) == [1, 2, 3, 4, 5]
    assert list(proteins[1].model_raw_data) == [2, 3, 4, 5, 6]

    def duplicate_outputs(args):
        write_outputs(args)
        with open(args[args.index('--out') + 1], 'a') as output:
            output.write('P1\t1.0\n')

    monkeypatch.setattr(Adhesin, 'run_command', duplicate_outputs)
    with pytest.raises(ValueError, match='duplicated'):
        Adhesin.extract_features([Protein('sp|P1|P1', 'AAAA'),
                                  Protein('sp|P2|P2', 'CCCC')], '', '',
                                 str(tmp_path), '')
