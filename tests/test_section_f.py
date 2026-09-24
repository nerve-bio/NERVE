import sys
import types

import numpy as np
import pytest

try:
    from code import Adhesin, Virulent_factor
except ModuleNotFoundError as error:
    if error.name != 'tensorflow':
        raise
    sys.modules['tensorflow'] = types.SimpleNamespace(keras=None)
    from code import Adhesin, Virulent_factor

from code.Checkpoint import compute_signature


class ProteinStub:
    def __init__(self, index):
        self.index = index
        self.p_ad = None
        self.p_vir = None

    def standardize(self, _means, _std_devs, _projection):
        return np.array([self.index], dtype=float)


@pytest.mark.parametrize('module,function,attribute', [
    (Adhesin, Adhesin.adhesin_predict, 'p_ad'),
    (Virulent_factor, Virulent_factor.virulent_factor_predict, 'p_vir'),
])
def test_inference_batch_keeps_protein_order_and_last_chunk(
        monkeypatch, module, function, attribute):
    calls = []

    class Model:
        def predict(self, data, batch_size, verbose):
            calls.append((len(data), batch_size, verbose))
            return data[:, 0] / 10000

    monkeypatch.setattr(module, 'keras', types.SimpleNamespace(
        models=types.SimpleNamespace(load_model=lambda _path: Model())))
    monkeypatch.setattr(module.np, 'load', lambda _path: np.array([1.0]))
    proteins = [ProteinStub(index) for index in range(2050)]
    assert function(proteins, '/unused', 64) is proteins
    assert calls == [(2048, 64, 0), (2, 64, 0)]
    assert [getattr(protein, attribute) for protein in proteins] == [
        index / 10000 for index in range(2050)]


def test_batch_one_remains_default(monkeypatch):
    observed = []

    class Model:
        def predict(self, data, batch_size, verbose):
            observed.append(batch_size)
            return np.zeros(len(data))

    monkeypatch.setattr(Adhesin, 'keras', types.SimpleNamespace(
        models=types.SimpleNamespace(load_model=lambda _path: Model())))
    monkeypatch.setattr(Adhesin.np, 'load', lambda _path: np.array([1.0]))
    Adhesin.adhesin_predict([ProteinStub(0)], '/unused')
    assert observed == [1]


def test_batch_setting_changes_checkpoint_signature():
    inputs = [{'name': 'input.fasta', 'sha256': 'abc'}]
    one = compute_signature([], {'inference_batch_size': 1}, input_info=inputs)
    larger = compute_signature([], {'inference_batch_size': 64}, input_info=inputs)
    assert one != larger
