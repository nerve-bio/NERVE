"""Short failure checks for the external boundaries of the pipeline."""

import io
import os
import subprocess
from types import SimpleNamespace

import pytest
import requests

from code import Function, Immunity, Quality_control, Subcellular, Utils


def test_external_command_keeps_paths_and_reports_failure(monkeypatch):
    arguments = ['tool', '--input', '/a path/input.fasta']
    observed = []
    def fake_run(argv, **kwargs):
        observed.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 7, '', 'bad input')
    monkeypatch.setattr(Utils.subprocess, 'run', fake_run)
    with pytest.raises(RuntimeError, match='tool exited 7: bad input'):
        Utils.run_command(arguments)
    assert observed[0][0] == arguments


def test_download_does_not_replace_valid_file_after_partial_failure(tmp_path, monkeypatch):
    target = tmp_path / 'proteome.fasta'
    target.write_text('>old\nAAAA\n')
    class PartialResponse:
        def raise_for_status(self):
            pass
        def iter_content(self, chunk_size):
            yield b'>new\nAAAA\n'
            raise requests.ConnectionError('stream interrupted')
        def close(self):
            pass
    monkeypatch.setattr(Quality_control.requests, 'get',
                        lambda *a, **k: PartialResponse())
    with pytest.raises(requests.ConnectionError, match='interrupted'):
        Quality_control.download_from_url_to_file('https://example.invalid', str(target))
    assert target.read_text() == '>old\nAAAA\n'
    assert not list(tmp_path.glob('.proteome-*'))


def test_download_rejects_empty_or_non_fasta(tmp_path, monkeypatch):
    class Response:
        def raise_for_status(self):
            pass
        def iter_content(self, chunk_size):
            yield b'not fasta\n'
        def close(self):
            pass
    monkeypatch.setattr(Quality_control.requests, 'get', lambda *a, **k: Response())
    with pytest.raises(ValueError, match='not FASTA'):
        Quality_control.download_from_url_to_file('https://example.invalid',
                                                   str(tmp_path / 'proteome.fasta'))
    assert not (tmp_path / 'proteome.fasta').exists()


def test_download_atomically_replaces_valid_fasta(tmp_path, monkeypatch):
    target = tmp_path / 'proteome.fasta'
    target.write_text('>old\nAAAA\n')
    class Response:
        def raise_for_status(self):
            pass
        def iter_content(self, chunk_size):
            yield b'>new\nACDE\n'
        def close(self):
            pass
    monkeypatch.setattr(Quality_control.requests, 'get', lambda *a, **k: Response())
    Quality_control.download_from_url_to_file('https://example.invalid', str(target))
    assert target.read_text() == '>new\nACDE\n'
    assert not list(tmp_path.glob('.proteome-*'))


def test_deepfri_failure_preserves_parent_directory_and_rejects_stale_file(tmp_path, monkeypatch):
    output = tmp_path / 'output path'
    output.mkdir()
    (output / '_MF_predictions.csv').write_text('stale')
    (output / '_MF_pred_scores.json').write_text('stale')
    original = os.getcwd()
    def fail(argv, *, cwd=None):
        assert cwd == str(tmp_path / 'deepfri')
        assert argv[-1] == str(output) + os.sep
        raise RuntimeError('DeepFRI exited 1')
    monkeypatch.setattr(Function, 'run_command', fail)
    with pytest.raises(RuntimeError, match='DeepFRI exited 1'):
        Function.deep_fri(str(tmp_path / 'input path.fasta'), str(tmp_path / 'deepfri'),
                          str(output))
    assert os.getcwd() == original
    assert not (output / '_MF_predictions.csv').exists()


def test_blast_rejects_missing_new_xml(tmp_path, monkeypatch):
    xml = tmp_path / 'old.xml'
    xml.write_text('stale')
    monkeypatch.setattr(Immunity, 'run_command', lambda *a, **k: ('', ''))
    with pytest.raises(RuntimeError, match='did not produce'):
        Immunity.run_blast('/input path.fasta', '/db path', 1e-10, str(xml), 1)
    assert not xml.exists()


def test_psortb_rejects_missing_duplicate_and_unknown_results():
    proteins = [SimpleNamespace(id='sp|P1|P1', accession='P1', localization=None),
                SimpleNamespace(id='sp|P2|P2', accession='P2', localization=None)]
    def entry(identifier):
        return (f'SeqID: {identifier}\nLocalization Scores:\n'
                'Cytoplasmic 8.96\nExtracellular 0.26\nFinal Prediction:\n')
    separator = '-------------------------------------------------------------------------------\n\n'
    with pytest.raises(ValueError, match='omitted 1 proteins'):
        Subcellular.parse_psortb_result(entry('P1'), proteins)
    with pytest.raises(ValueError, match='duplicate protein'):
        Subcellular.parse_psortb_result(separator.join([entry('P1'), entry('P1')]), proteins)
    with pytest.raises(ValueError, match='unknown protein'):
        Subcellular.parse_psortb_result(entry('P3'), proteins)
    result = Subcellular.parse_psortb_result(separator.join([entry('P1'), entry('P2')]), proteins)
    assert result['P1'][0].localization == 'Cytoplasmic'


def test_psortb_service_failure_leaves_proteins_unmodified(tmp_path, monkeypatch):
    input_path = tmp_path / 'input.fasta'
    input_path.write_text('>sp|P1|P1\nAAAA\n')
    protein = SimpleNamespace(id='sp|P1|P1', accession='P1', localization=None)
    class Response:
        def __enter__(self):
            return io.BytesIO(b'{"error":"failed","result":""}')
        def __exit__(self, *args):
            pass
    monkeypatch.setattr(Subcellular.urllib.request, 'urlopen', lambda *a, **k: Response())
    with pytest.raises(RuntimeError, match='service error'):
        Subcellular.psortb([protein], str(tmp_path), 'n', str(input_path))
    assert protein.localization is None
