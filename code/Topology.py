#!/usr/local/bin/python
"""Runs prediction on protein topology with TMHMM"""

import logging
import multiprocessing as mp
import os
import signal

import tmhmm
from tmhmm.api import DEFAULT_MODEL as _DEFAULT_TMHMM_MODEL
from tmhmm.model import parse as _parse_tmhmm_model
from code.Utils import Progress


_MODEL = None


def _tmhmm_model():
    global _MODEL
    if _MODEL is None:
        _, _MODEL = _parse_tmhmm_model(_DEFAULT_TMHMM_MODEL)
    return _MODEL


def _worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    _tmhmm_model()


def _predict_path(sequence):
    return tmhmm.predict(sequence, _tmhmm_model(), compute_posterior=False)


def _count_helices(annotation):
    return sum(annotation[i] in ('i', 'o') and annotation[i + 1] == 'M'
               for i in range(len(annotation) - 1))


def tmhelices(list_of_proteins, working_dir, threads=1)->list:
    "Runs Tmhelices"
    sequences = [protein.sequence for protein in list_of_proteins]
    workers = min(max(1, int(threads)), max(1, len(sequences)))
    progress = Progress('topology', len(sequences))
    if workers > 1:
        chunksize = max(1, min(64, len(sequences) // (workers * 8) or 1))
        pool = mp.get_context('fork').Pool(workers, initializer=_worker_init)
        try:
            annotations = []
            for annotation in pool.imap(_predict_path, sequences, chunksize=chunksize):
                annotations.append(annotation)
                progress.update(len(annotations))
        finally:
            pool.terminate()
            pool.join()
    else:
        _tmhmm_model()
        annotations = []
        for sequence in sequences:
            annotations.append(_predict_path(sequence))
            progress.update(len(annotations))

    for protein, annotation in zip(list_of_proteins, annotations):
        protein.tmhmm_seq = annotation
        protein.transmembrane_doms = _count_helices(annotation)
    return list_of_proteins
