#!/usr/bin/env python3
"""Run real epitopepredict objects through serial and auto ProcessPool paths."""
import os, shutil, tempfile
from code.Protein import Protein
from code.Epitope import epitope, available_cpu_count

def candidates():
    sequences = ['MKTIIALSYIFCLVFADYKDDDDK', 'GILGFVFTLTVPSERKTYQGSYG',
                 'LLFGYPVYVLVTRGQYGRSLLQG', 'MTEYKLVVVGAGGVGKSALTIQL']
    result = []
    for index, sequence in enumerate(sequences):
        protein = Protein(f'sp|P{index}|P{index} probe', sequence)
        protein.score = 1.0
        result.append(protein)
    return result

def signature(proteins):
    return [(p.MHC1_binders, p.MHC2_binders, p.MHC1_pb_binders, p.MHC2_pb_binders)
            for p in proteins]

def main():
    root = tempfile.mkdtemp(prefix='nerve-epitope-check-')
    try:
        serial = epitope(candidates(), os.path.join(root, 'serial'), 9, 11, 1, 1,
                         0.0, 'False', 3, 1)
        parallel = epitope(candidates(), os.path.join(root, 'parallel'), 9, 11, 1, 1,
                           0.0, 'False', 3, 0)
        if signature(serial) != signature(parallel):
            raise SystemExit('FAIL: serial and parallel epitope attributes differ')
        print(f'PASS: real predictor output equivalent; visible_cpus={available_cpu_count()}')
    finally:
        shutil.rmtree(root)

if __name__ == '__main__': main()
