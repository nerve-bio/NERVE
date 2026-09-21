"""Failure cases for the comparison tools used as validation gates."""

import pickle
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from harness.compare_states import compare_snapshots


ROOT = Path(__file__).resolve().parents[1]


def test_duplicate_snapshot_keys_fail(tmp_path):
    path = tmp_path / 'state.pkl'
    proteins = [SimpleNamespace(accession='P1', id='P1', score=1),
                SimpleNamespace(accession='P1', id='P1', score=2)]
    path.write_bytes(pickle.dumps(proteins))
    with pytest.raises(ValueError, match='duplicate protein key'):
        compare_snapshots(path, path)


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


def test_missing_epitope_outputs_fail(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir()
    b.mkdir()
    result = subprocess.run(
        [sys.executable, str(ROOT / 'harness/validation/compare_epitope_outputs.py'),
         str(a), str(b)], text=True, capture_output=True,
    )
    assert result.returncode != 0
