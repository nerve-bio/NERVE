#!/usr/bin/python3
"""
    Equivalence and unit tests for the identity and mhcpep-matching helpers.

    These run without Docker and without the external predictors: they exercise
    pure-Python helpers against the real mhcpep database. Run with:
        pytest tests/test_equivalence.py

    The mhcpep test is the important one. Replacing the exhaustive epitope scan
    with a k-mer index is only legitimate if it produces the identical list, in
    the identical order, so that is asserted directly against the naive
    implementation rather than argued for in a comment.
"""

import os
import sys
import random

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Biopython dropped the deprecated command-line wrappers in 1.80; the container
# pins 1.78, where they still exist. Nothing here invokes blastp -- the import is
# only needed for code.Immunity to load -- so stub it when running against a
# newer Biopython outside the container.
try:
    from Bio.Blast.Applications import NcbiblastpCommandline  # noqa: F401
except ImportError:
    import types
    _stub = types.ModuleType('Bio.Blast.Applications')
    _stub.NcbiblastpCommandline = object
    sys.modules['Bio.Blast.Applications'] = _stub

from code.Protein import Protein
from code.Immunity import build_kmer_index, compare_with_mhcpep
from code.Utils import (derive_accession, normalize_header, build_id_index,
                        match_protein)

NERVE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MHCPEP_SAPIENS = os.path.join(NERVE_DIR, 'database/mhcpep/mhcpep_sapiens.csv')

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def naive_mhcpep(list_of_proteins, peptides, shared_attr, out_attr):
    """The original exhaustive scan, kept here as the reference implementation."""
    for protein in list_of_proteins:
        for seq in getattr(protein, shared_attr):
            for peptide in peptides:
                getattr(protein, out_attr).extend(Protein.peptide_comparison(seq, peptide))


def make_shared_peptide(query, start_pos, length):
    """Build a shared-peptide record shaped like hsp_match_parser's output."""
    return {'query': query, 'match': query[start_pos:start_pos + length],
            'start_pos': start_pos}


def build_case(peptides, seed):
    """A protein whose shared peptides are a mix of real epitope substrings
    (which must match) and random sequence (which mostly must not)."""
    rng = random.Random(seed)
    real = [p for p in peptides if isinstance(p, str) and len(p) >= 12]
    records = []
    for _ in range(15):
        source = rng.choice(real)
        start = rng.randrange(0, max(1, len(source) - 9))
        # embed a real epitope window inside a longer synthetic query
        flank_l = ''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randrange(0, 12)))
        flank_r = ''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randrange(0, 12)))
        query = flank_l + source[start:start + 9] + flank_r
        records.append(make_shared_peptide(query, len(flank_l), 9))
    for _ in range(10):
        query = ''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randrange(12, 30)))
        records.append(make_shared_peptide(query, rng.randrange(0, 3), 9))
    return records


def test_mhcpep_index_is_exactly_equivalent():
    """The k-mer index must reproduce the exhaustive scan exactly, in order."""
    peptides = list(pd.read_csv(MHCPEP_SAPIENS, skipinitialspace=True)['Epitope.2'])

    for seed in range(5):
        shared = build_case(peptides, seed)

        reference = Protein('sp|REF|REF reference', 'A' * 100)
        reference.list_of_shared_human_peps = shared
        naive_mhcpep([reference], peptides, 'list_of_shared_human_peps',
                     'list_of_peptides_from_comparison_with_mhcpep_sapiens')

        indexed = Protein('sp|IDX|IDX indexed', 'A' * 100)
        indexed.list_of_shared_human_peps = shared
        compare_with_mhcpep([indexed], peptides, 'list_of_shared_human_peps',
                            'list_of_peptides_from_comparison_with_mhcpep_sapiens')

        expected = reference.list_of_peptides_from_comparison_with_mhcpep_sapiens
        actual = indexed.list_of_peptides_from_comparison_with_mhcpep_sapiens
        assert actual == expected, f'seed {seed}: {len(actual)} matches vs {len(expected)} expected'
        # a test that matched nothing on either side would pass vacuously
        assert len(expected) > 0, f'seed {seed} produced no matches; the case is not exercising anything'


def test_kmer_index_skips_only_non_matching_peptides():
    """Every peptide containing the window must be indexed under it, in file order."""
    peptides = ['AAACCCGGG', 'CCCGGGTTT', 'YYYYYYYYY', None, float('nan'), 'CCC']
    index = build_kmer_index(peptides, 6)
    assert index['CCCGGG'] == ['AAACCCGGG', 'CCCGGGTTT']
    assert 'YYYYYY' in index
    # too short to contain a 6-mer, and the non-string entries, are skipped
    assert all('CCC' != p for p in index.get('CCCGGG', []))


def test_derive_accession():
    assert derive_accession('sp|P0C277|FOLD_NEIMB Bifunctional protein') == 'P0C277'
    assert derive_accession('tr|Q0P9Y0|Q0P9Y0_CAMJE Putative') == 'Q0P9Y0'
    # NCBI style with more than three fields used to yield None
    assert derive_accession('gi|12345|ref|NP_001|') == '12345'
    # bare in-house identifiers fall back to the first token
    assert derive_accession('LBRM2904_000010100 hypothetical protein') == 'LBRM2904_000010100'
    assert derive_accession('single') == 'single'


def test_normalize_header_is_unique_and_safe():
    seen = set()
    first, acc1 = normalize_header('sp|P0C277|FOLD_NEIMB Bifunctional protein', seen)
    assert first == 'sp|P0C277|P0C277 Bifunctional protein'
    assert acc1 == 'P0C277'
    # a colliding accession must be disambiguated, not silently merged
    _, acc2 = normalize_header('sp|P0C277|OTHER Another protein', seen)
    assert acc2 == 'P0C277_2'
    # path separators cannot survive: the accession names an output directory
    _, acc3 = normalize_header('weird/id|with|slashes description', seen)
    assert '/' not in acc3


def test_exact_matching_beats_substring_ambiguity():
    """A short identifier must not bind to the longer protein containing it."""
    short = Protein('sp|P1|P1 short', 'AAAA')
    longer = Protein('sp|P10|P10 longer', 'CCCC')
    index = build_id_index([longer, short])
    assert match_protein('sp|P1|P1', index, [longer, short]) is short
    assert match_protein('P1', index, [longer, short]) is short
    assert match_protein('P10', index, [longer, short]) is longer


def test_unresolvable_identifier_returns_none():
    p = Protein('sp|P1|P1 only', 'AAAA')
    index = build_id_index([p])
    assert match_protein('sp|NOPE|NOPE', index, None) is None


def test_accession_is_never_none():
    for header in ('bare_identifier', 'gi|1|ref|NP_1|', 'sp|P0C277|X d', 'a b c'):
        assert Protein(header, 'AAAA').accession
