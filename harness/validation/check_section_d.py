#!/usr/bin/env python3
"""Check Section D run records and terminal outputs."""

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ('vaccine_candidates.csv', 'discarded_proteins.csv',
           'selection_reasons.csv')
REFERENCE_HASHES = {
    'vaccine_candidates.csv': '9ef1defbe12dafe8161300b212f17ea210dc2b963370577bb2107b6e08aba5f0',
    'discarded_proteins.csv': 'c0a01d68caf8589e4eed8e20aaf1cdea1f3c51b7510820e741e65c2b3770bb09',
}


def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def host_path(path):
    prefix = '/workdir/'
    if not path.startswith(prefix):
        raise ValueError(f'Unexpected container path: {path}')
    return ROOT / path[len(prefix):]


def check_run(directory, commit, attempts, last_stage):
    manifest = json.loads((directory / 'run_manifest.json').read_text())
    assert manifest['status'] == 'success', directory
    assert manifest['source']['commit'] == commit, directory
    assert manifest['run_signature'], directory
    assert manifest['attempt'] == attempts, directory
    assert len(manifest['attempts']) == attempts, directory
    assert manifest['attempts'][-1]['status'] == 'success', directory
    if attempts == 2:
        assert manifest['attempts'][0]['status'] == 'interrupted', directory
    assert manifest['last_completed_stage'] == last_stage, directory
    assert manifest['attempt_wall_seconds'] >= 0, directory
    assert {stage['ordinal'] for stage in manifest['stages']} == (
        {1, 2, 3, 4, 6, 11} if last_stage == 11 else set(range(1, 13))), directory
    assert manifest['qc']['primary']['kept'] > 0, directory
    assert len(manifest['inputs']) == (1 if last_stage == 11 else 2), directory
    for item in manifest['inputs']:
        path = host_path(item['path'])
        assert path.stat().st_size == item['bytes'], path
        assert digest(path) == item['sha256'], path
    for name in OUTPUTS:
        path = directory / name
        record = manifest['outputs'][name]
        assert host_path(record['path']) == path, path
        assert digest(path) == record['sha256'], path
    with open(directory / 'selection_reasons.csv', newline='') as handle:
        reasons = list(csv.DictReader(handle))
    with open(directory / 'vaccine_candidates.csv', newline='') as handle:
        selected = {row['uniprot_accession_code'] for row in csv.DictReader(handle)}
    with open(directory / 'discarded_proteins.csv', newline='') as handle:
        discarded = {row['uniprot_accession_code'] for row in csv.DictReader(handle)}
    assert selected.isdisjoint(discarded), directory
    assert {row['accession'] for row in reasons if row['selected'] == 'True'} == selected, directory
    assert {row['accession'] for row in reasons if row['selected'] == 'False'} == discarded, directory
    assert dict(Counter(row['first_rule'] for row in reasons)) == (
        manifest['selection_reason_counts']), directory
    assert manifest['candidate_count'] == len(selected), directory
    assert manifest['discarded_count'] == len(discarded), directory
    log = (directory / 'run.log').read_text()
    assert 'Stage 1 quality_control:' in log, directory
    assert 'End NERVE computation successfully' in log, directory
    assert '10% done' not in log, directory
    assert not (directory / 'debug.log').exists(), directory
    assert not list(directory.glob('.manifest-*')), directory
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-root', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--image', required=True)
    args = parser.parse_args()
    image = json.loads(subprocess.check_output(
        ['docker', 'image', 'inspect', args.image], text=True))[0]
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert image['Config']['Labels']['org.opencontainers.image.revision'] == commit
    runs = {
        'checkpoint/resume': check_run(args.run_root / 'checkpoint/resume', commit, 2, 11),
        'checkpoint/fresh': check_run(args.run_root / 'checkpoint/fresh', commit, 1, 11),
        'cjejuni/reference-t1-ew1': check_run(
            args.run_root / 'cjejuni/reference-t1-ew1', commit, 1, 12),
    }
    mismatch = args.run_root / 'checkpoint/mismatch'
    rejected = json.loads((mismatch / 'run_manifest.json').read_text())
    assert rejected['status'] in ('running', 'interrupted'), mismatch
    assert rejected['attempt'] == 1, mismatch
    assert not (mismatch / 'vaccine_candidates.csv').exists(), mismatch
    assert 'Checkpoint does not match' in (mismatch / 'expected-rejection.log').read_text(), mismatch
    current = args.run_root / 'cjejuni/reference-t1-ew1'
    for name, expected in REFERENCE_HASHES.items():
        assert digest(args.reference / name) == expected, name
        assert digest(current / name) == expected, name
    report = {
        'status': 'passed', 'commit': commit, 'image_id': image['Id'],
        'reference': str(args.reference.resolve()),
        'runs': {name: {'status': run['status'], 'attempts': run['attempt'],
                        'candidate_count': run['candidate_count']}
                 for name, run in runs.items()},
        'reference_output_sha256': REFERENCE_HASHES,
    }
    (args.run_root / 'section_d_gate.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Section D manifest, logging, selection, and predecessor checks passed.')


if __name__ == '__main__':
    main()
