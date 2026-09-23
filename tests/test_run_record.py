"""Run records and selection explanations do not alter scientific output."""

import json
from pathlib import Path
from types import SimpleNamespace

from code.Run_record import RunRecord, input_records
from code.Select import select, write_selection_reasons


def test_atomic_run_record_tracks_stage_failure_and_resume(tmp_path):
    fasta = tmp_path / 'input.fasta'
    fasta.write_text('>P1\nACDE\n')
    args = SimpleNamespace(_asdict=lambda: {'gram': 'n', 'debug': 'False'})
    inputs = input_records([fasta])
    record = RunRecord(tmp_path, args, inputs, {'commit': 'abc'}, 'signature')
    record.stage(1, 'quality_control', 0.25, 1)
    record.finish('interrupted', 1.5, error='interrupted')
    saved = json.loads((tmp_path / 'run_manifest.json').read_text())
    assert saved['status'] == 'interrupted'
    assert saved['last_completed_stage'] == 1
    assert saved['inputs'][0]['sha256'] == inputs[0]['sha256']
    resumed = RunRecord(tmp_path, args, inputs, {'commit': 'abc'}, 'signature', resume=True)
    resumed.restored_through(1)
    resumed.finish('success', 2.0, candidate_count=1)
    saved = json.loads((tmp_path / 'run_manifest.json').read_text())
    assert saved['attempt'] == 2 and saved['status'] == 'success'
    assert not list(tmp_path.glob('.manifest-*'))


def test_selection_sidecar_agrees_with_existing_decisions(tmp_path):
    def protein(accession, location, adhesion):
        return SimpleNamespace(
            accession=accession,
            localization=[SimpleNamespace(localization=location, reliability=9.0)],
            p_ad=adhesion, p_vir=None, transmembrane_doms=0,
            original_sequence_if_razor=None, sapiens_peptides_sum=0.0,
            list_of_peptides_from_comparison_with_mhcpep_sapiens=[],
            mouse_peptides_sum=0.0, list_of_peptides_from_comparison_with_mhcpep_mouse=[])
    proteins = [protein('P1', 'OuterMembrane', .9),
                protein('P2', 'Cytoplasmic', .9),
                protein('P3', 'OuterMembrane', .1)]
    args = SimpleNamespace(select='True', transmemb_doms_limit=3,
                           padlimit=.5, mouse='False', mouse_peptides_sum_limit=.15,
                           virlimit=.5, virulent='False', razor='False')
    selected = select(proteins, 3, .5, 'False', .15, .5, 'False', 'False')
    counts = write_selection_reasons(proteins, args, tmp_path / 'selection_reasons.csv')
    rows = (tmp_path / 'selection_reasons.csv').read_text()
    assert [p.accession for p in selected] == ['P1']
    assert counts == {'selected': 1, 'cytoplasmic': 1, 'adhesin_probability': 1}
    assert 'P1,True,selected' in rows and 'P2,False,cytoplasmic' in rows
    args.select = 'False'
    assert write_selection_reasons(proteins, args, tmp_path / 'disabled.csv') == {
        'selection_disabled': 3}
