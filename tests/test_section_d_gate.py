"""Section D gate checks saved output, not just process exit."""

import json

import pytest

from harness.validation import check_section_d


def test_run_record_gate_checks_selection_and_hashes(tmp_path, monkeypatch):
    monkeypatch.setattr(check_section_d, 'ROOT', tmp_path)
    directory = tmp_path / 'run'
    directory.mkdir()
    fasta = tmp_path / 'input.fasta'
    fasta.write_text('>P1\nACDE\n>P2\nACDE\n')
    (directory / 'vaccine_candidates.csv').write_text('uniprot_accession_code\nP1\n')
    (directory / 'discarded_proteins.csv').write_text('uniprot_accession_code\nP2\n')
    (directory / 'selection_reasons.csv').write_text(
        'accession,selected,first_rule,observed,cutoff\n'
        'P1,True,selected,,\nP2,False,cytoplasmic,Cytoplasmic,non-cytoplasmic\n')
    (directory / 'run.log').write_text(
        'Stage 1 quality_control: starting\nEnd NERVE computation successfully\n')
    outputs = {name: {'path': '/workdir/run/' + name,
                      'sha256': check_section_d.digest(directory / name)}
               for name in check_section_d.OUTPUTS}
    manifest = {
        'status': 'success', 'source': {'commit': 'commit'}, 'run_signature': 'signature',
        'attempt': 1, 'attempts': [{'status': 'success'}], 'attempt_wall_seconds': 1.0,
        'last_completed_stage': 11,
        'stages': [{'ordinal': number} for number in (1, 2, 3, 4, 6, 11)],
        'qc': {'primary': {'kept': 2}},
        'inputs': [{'path': '/workdir/input.fasta', 'bytes': fasta.stat().st_size,
                    'sha256': check_section_d.digest(fasta)}],
        'outputs': outputs, 'candidate_count': 1, 'discarded_count': 1,
        'selection_reason_counts': {'selected': 1, 'cytoplasmic': 1},
    }
    (directory / 'run_manifest.json').write_text(json.dumps(manifest))
    assert check_section_d.check_run(directory, 'commit', 1, 11) == manifest
    (directory / 'selection_reasons.csv').write_text(
        'accession,selected,first_rule,observed,cutoff\nP1,False,cytoplasmic,,\n')
    manifest['outputs']['selection_reasons.csv']['sha256'] = check_section_d.digest(
        directory / 'selection_reasons.csv')
    (directory / 'run_manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(AssertionError):
        check_section_d.check_run(directory, 'commit', 1, 11)
