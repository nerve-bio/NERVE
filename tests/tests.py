"""Docker integration checks for the optimized pipeline.

Run after building the optimized image and starting PSORTb on nerve-network:
    NERVE_TEST_IMAGE=nerve:v0.1.0 pytest -q tests/tests.py -p no:debugging
Each test mounts a fresh output directory and leaves the repository untouched.
"""

import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
IMAGE = os.environ.get('NERVE_TEST_IMAGE', 'nerve:v0.1.0')


def run_image(output, *args):
    command = [
        'docker', 'run', '--rm', '--network', 'nerve-network',
        '-v', f'{ROOT / "tests" / "input_data"}:/inputs:ro',
        '-v', f'{output}:/outputs', IMAGE, *args,
    ]
    return subprocess.run(command, text=True, capture_output=True)


def test_help(tmp_path):
    result = run_image(tmp_path, '--help')
    assert result.returncode == 0, result.stderr
    assert '--checkpoint' in result.stdout


@pytest.mark.parametrize('extra', [
    ('--select', 'False'),
    ('--select', 'True'),
    ('--mouse', 'True', '--select', 'True'),
])
def test_small_pipeline(tmp_path, extra):
    result = run_image(
        tmp_path, '--proteome1', '/inputs/proteome1.fasta', '--gram', 'n',
        '--working_dir', '/outputs', '--annotation', 'False',
        '--razor', 'False', '--virulent', 'False', '--epitopes', 'False',
        '--threads', '1', *extra,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'End NERVE computation successfully.' in result.stdout
    candidates = tmp_path / 'vaccine_candidates.csv'
    discarded = tmp_path / 'discarded_proteins.csv'
    assert candidates.is_file() and discarded.is_file()
    for path in (candidates, discarded):
        frame = pd.read_csv(path)
        assert {'id', 'uniprot_accession_code', 'score'} <= set(frame.columns)


def test_epitope_artifacts(tmp_path):
    result = run_image(
        tmp_path, '--proteome1', '/inputs/proteome5.fasta', '--gram', 'n',
        '--working_dir', '/outputs', '--annotation', 'False',
        '--razor', 'False', '--mouse', 'False', '--virulent', 'False',
        '--select', 'True', '--epitopes', 'True', '--ep_plots', 'False',
        '--epitope_workers', '1',
    )
    assert result.returncode == 0, result.stdout + result.stderr
    epitope = tmp_path / 'epitope'
    assert epitope.is_dir()
    assert list(epitope.rglob('mhci_epitopes_*.csv'))
    assert list(epitope.rglob('mhcii_epitopes_*.csv'))
