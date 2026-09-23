"""Run manifest."""

import hashlib
import json
import os
import platform
import tempfile
from datetime import datetime, timezone


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def input_records(paths):
    return [{'path': os.path.abspath(path), 'name': os.path.basename(path),
             'bytes': os.path.getsize(path), 'sha256': file_hash(path)}
            for path in paths if path]


def atomic_json(path, data):
    directory = os.path.dirname(path)
    with tempfile.NamedTemporaryFile('w', dir=directory, prefix='.manifest-',
                                     delete=False) as output:
        temporary = output.name
        try:
            json.dump(data, output, indent=2, sort_keys=True)
            output.write('\n')
            output.flush()
            os.fsync(output.fileno())
            os.chmod(temporary, 0o644)
        except BaseException:
            os.unlink(temporary)
            raise
    os.replace(temporary, path)


class RunRecord:
    def __init__(self, directory, args, inputs, source, signature, resume=False):
        self.path = os.path.join(directory, 'run_manifest.json')
        if resume:
            with open(self.path) as handle:
                self.data = json.load(handle)
            if self.data['status'] == 'running':
                self.data.setdefault('attempts', []).append({
                    'number': self.data['attempt'], 'status': 'interrupted',
                    'started_utc': self.data['attempt_started_utc'],
                    'ended_utc': None, 'wall_seconds': None,
                    'error': 'Previous process stopped without a final status'})
            self.data['attempt'] += 1
        else:
            self.data = {
                'schema': 1, 'started_utc': utc_now(), 'attempt': 1,
                'source': source, 'inputs': inputs,
                'host': {'architecture': platform.machine(),
                         'logical_cpus': os.cpu_count()},
                'image': {'id': os.environ.get('NERVE_IMAGE_ID'),
                          'base_digest': os.environ.get('NERVE_BASE_IMAGE_DIGEST'),
                          'psortb_id': os.environ.get('NERVE_PSORTB_IMAGE_ID')},
                'cpu_limit': os.environ.get('NERVE_CPU_LIMIT'),
                'memory_limit_bytes': os.environ.get('NERVE_MEMORY_LIMIT_BYTES'),
                'stages': [], 'attempts': [],
            }
        self.data.update({'arguments': args._asdict(), 'run_signature': signature,
                          'attempt_started_utc': utc_now(), 'attempt_ended_utc': None,
                          'attempt_wall_seconds': None, 'status': 'running', 'error': None})
        self.save()

    def save(self):
        atomic_json(self.path, self.data)

    def stage(self, ordinal, name, seconds, count):
        self.data['stages'] = [stage for stage in self.data['stages']
                               if stage['ordinal'] != ordinal]
        self.data['stages'].append({'ordinal': ordinal, 'name': name,
                                    'seconds': round(seconds, 4), 'proteins': count})
        self.data['last_completed_stage'] = ordinal
        self.save()

    def restored_through(self, ordinal):
        self.data['last_completed_stage'] = max(
            ordinal, self.data.get('last_completed_stage', 0))
        self.save()

    def finish(self, status, elapsed, error=None, **fields):
        self.data.update(fields)
        self.data.update({'status': status, 'error': error,
                          'attempt_ended_utc': utc_now(),
                          'attempt_wall_seconds': round(elapsed, 4)})
        self.data['attempts'].append({
            'number': self.data['attempt'], 'status': status,
            'started_utc': self.data['attempt_started_utc'],
            'ended_utc': self.data['attempt_ended_utc'],
            'wall_seconds': self.data['attempt_wall_seconds'], 'error': error})
        self.save()
