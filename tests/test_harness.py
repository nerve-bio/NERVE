"""Failure cases for the comparison tools used as validation gates."""

import pickle
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from harness.compare_states import compare_snapshots
from code.Checkpoint import Checkpointer, compute_signature
from code.Subcellular import Localization


ROOT = Path(__file__).resolve().parents[1]


def test_duplicate_snapshot_keys_fail(tmp_path):
    path = tmp_path / 'state.pkl'
    proteins = [SimpleNamespace(accession='P1', id='P1', score=1),
                SimpleNamespace(accession='P1', id='P1', score=2)]
    path.write_bytes(pickle.dumps(proteins))
    with pytest.raises(ValueError, match='duplicate protein key'):
        compare_snapshots(path, path)


def test_missing_and_changed_snapshot_proteins_fail(tmp_path):
    a, b = tmp_path / 'a.pkl', tmp_path / 'b.pkl'
    a.write_bytes(pickle.dumps([SimpleNamespace(accession='P1', score=1)]))
    b.write_bytes(pickle.dumps([SimpleNamespace(accession='P2', score=1)]))
    result = compare_snapshots(a, b)
    assert result['only_a'] == ['P1'] and result['only_b'] == ['P2']
    b.write_bytes(pickle.dumps([SimpleNamespace(accession='P1', score=2)]))
    assert compare_snapshots(a, b)['per_attr']['score']['count'] == 1


def test_only_secondary_localization_score_ties_are_canonical(tmp_path):
    a, b = tmp_path / 'a.pkl', tmp_path / 'b.pkl'
    def save(path, names):
        locations = [Localization(name, score) for name, score in names]
        path.write_bytes(pickle.dumps([SimpleNamespace(accession='P1', localization=locations)]))
    first = [('Cytoplasmic', 8.96), ('Extracellular', 0.26),
             ('Periplasmic', 0.26), ('OuterMembrane', 0.01)]
    save(a, first)
    save(b, [first[0], first[2], first[1], first[3]])
    assert compare_snapshots(a, b)['per_attr']['localization']['count'] == 1
    assert not compare_snapshots(a, b, canonical_localization_ties=True)['per_attr']
    save(b, [('Extracellular', 8.96), *first[1:]])
    assert compare_snapshots(a, b, canonical_localization_ties=True)['per_attr']
    save(b, [first[0], ('Periplasmic', 0.27), first[1], first[3]])
    assert compare_snapshots(a, b, canonical_localization_ties=True)['per_attr']


def compare_candidates(a, b):
    return subprocess.run(
        [sys.executable, str(ROOT / 'harness/compare_candidates.py'), str(a), str(b)],
        text=True, capture_output=True,
    )


def test_candidate_columns_can_differ(tmp_path):
    a, b = tmp_path / 'a.csv', tmp_path / 'b.csv'
    pd.DataFrame({'uniprot_accession_code': ['P1'], 'score': [1]}).to_csv(a, index=False)
    pd.DataFrame({'id': ['P1'], 'score': [1]}).to_csv(b, index=False)
    assert compare_candidates(a, b).returncode == 0
    pd.DataFrame({'id': ['P2'], 'score': [1]}).to_csv(b, index=False)
    assert compare_candidates(a, b).returncode == 1


def test_duplicate_and_missing_candidate_ids_fail(tmp_path):
    a, b = tmp_path / 'a.csv', tmp_path / 'b.csv'
    pd.DataFrame({'id': ['P1'], 'score': [1]}).to_csv(a, index=False)
    pd.DataFrame({'id': ['P1', 'P1'], 'score': [1, 2]}).to_csv(b, index=False)
    assert compare_candidates(a, b).returncode != 0
    pd.DataFrame({'id': [None], 'score': [1]}).to_csv(b, index=False)
    assert compare_candidates(a, b).returncode != 0


def test_changed_missing_and_duplicate_candidates_explain_failure(tmp_path):
    a, b = tmp_path / 'a.csv', tmp_path / 'b.csv'
    pd.DataFrame({'id': ['P1'], 'score': [1]}).to_csv(a, index=False)
    for ids, scores, reason in [(['P1'], [2], 'DIFFERENT'),
                                (['P2'], [1], 'ONLY in'),
                                (['P1', 'P1'], [1, 1], 'duplicate candidate identifier')]:
        pd.DataFrame({'id': ids, 'score': scores}).to_csv(b, index=False)
        result = compare_candidates(a, b)
        assert result.returncode != 0
        assert reason in result.stdout + result.stderr


def test_checkpoint_rejects_wrong_input_parameters_and_source(tmp_path):
    source = tmp_path / 'input.fasta'
    source.write_text('>P1\nAAAA\n')
    signature = compute_signature([source], {'gram': 'n'}, {'tree': 'one'})
    Checkpointer(tmp_path, signature, enabled=True).save(1, 'qc', [])
    alternatives = []
    source.write_text('>P1\nAAAT\n')
    alternatives.append(compute_signature([source], {'gram': 'n'}, {'tree': 'one'}))
    source.write_text('>P1\nAAAA\n')
    alternatives.append(compute_signature([source], {'gram': 'p'}, {'tree': 'one'}))
    alternatives.append(compute_signature([source], {'gram': 'n'}, {'tree': 'two'}))
    for wrong in alternatives:
        with pytest.raises(ValueError, match='Checkpoint does not match'):
            Checkpointer(tmp_path, wrong, resume=True)


def test_missing_epitope_outputs_fail(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir()
    b.mkdir()
    result = subprocess.run(
        [sys.executable, str(ROOT / 'harness/validation/compare_epitope_outputs.py'),
         str(a), str(b)], text=True, capture_output=True,
    )
    assert result.returncode != 0


def test_changed_and_missing_epitope_files_fail(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir()
    b.mkdir()
    pa, pb = a / 'epitopes.csv', b / 'epitopes.csv'
    pd.DataFrame({'allele': ['A'], 'peptide': ['AAAA'], 'score': [1]}).to_csv(pa, index=False)
    command = [sys.executable, str(ROOT / 'harness/validation/compare_epitope_outputs.py'),
               str(a), str(b)]
    missing = subprocess.run(command, text=True, capture_output=True)
    assert missing.returncode != 0
    assert 'must contain CSV files' in missing.stderr
    pd.DataFrame({'allele': ['A'], 'peptide': ['AAAA'], 'score': [2]}).to_csv(pb, index=False)
    changed = subprocess.run(command, text=True, capture_output=True)
    assert changed.returncode != 0 and 'DIFFERENT' in changed.stdout
    pd.DataFrame({'allele': ['B'], 'peptide': ['BBBB'], 'score': [1]}).to_csv(
        a / 'additional.csv', index=False)
    extra = subprocess.run(command, text=True, capture_output=True)
    assert extra.returncode != 0 and 'Only in A' in extra.stdout
