"""Crash-safe, input- and parameter-keyed pipeline checkpoints."""

import hashlib
import json
import logging
import os
import pickle


SCHEMA_VERSION = 1


def _file_digest(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def compute_signature(input_paths, parameters, source_identity=None):
    """Return a stable SHA-256 key for inputs, effective parameters and source."""
    inputs = []
    for path in input_paths:
        if path:
            inputs.append({'name': os.path.basename(path), 'sha256': _file_digest(path)})
    payload = {
        'schema': SCHEMA_VERSION,
        'inputs': inputs,
        'parameters': parameters,
        'source': source_identity or {},
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'), default=str).encode()
    return hashlib.sha256(encoded).hexdigest()


class Checkpointer:
    FILENAME = 'checkpoint.pkl'

    def __init__(self, working_dir, signature, enabled=False, resume=False):
        self.enabled = bool(enabled or resume)
        self.signature = signature
        self.path = os.path.join(working_dir, self.FILENAME)
        self.bundle = None
        if resume:
            self.bundle = self._load()

    def _load(self):
        if not os.path.isfile(self.path):
            logging.info('[checkpoint] no checkpoint found; starting from stage 1')
            return None
        with open(self.path, 'rb') as handle:
            bundle = pickle.load(handle)
        if bundle.get('schema') != SCHEMA_VERSION:
            raise ValueError('Checkpoint schema is incompatible with this NERVE version')
        if bundle.get('signature') != self.signature:
            raise ValueError('Checkpoint does not match the input, parameters, or source tree')
        logging.info('[checkpoint] resuming after stage %s (%s)',
                     bundle['ordinal'], bundle['name'])
        return bundle

    def done(self, ordinal):
        return self.bundle is not None and ordinal <= self.bundle['ordinal']

    def restored_state(self):
        return self.bundle.get('state') if self.bundle else None

    def save(self, ordinal, name, state):
        if not self.enabled:
            return
        bundle = {'schema': SCHEMA_VERSION, 'signature': self.signature,
                  'ordinal': ordinal, 'name': name, 'state': state}
        temporary = self.path + '.tmp'
        with open(temporary, 'wb') as handle:
            pickle.dump(bundle, handle, protocol=pickle.HIGHEST_PROTOCOL)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, self.path)
        self.bundle = bundle

    def finish(self):
        if self.enabled and os.path.isfile(self.path):
            os.remove(self.path)
