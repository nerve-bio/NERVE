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
from types import SimpleNamespace

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
from code.Immunity import build_kmer_index, compare_with_mhcpep, ordered_hsps
from code.Quality_control import quality_control, qc_sequence
from code.Checkpoint import Checkpointer, compute_signature
from code.Utils import (derive_accession, normalize_header, build_id_index,
                        match_protein, epitope_output_identifier)

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


def test_ambiguous_substring_fallback_fails_closed():
    first = Protein('legacy-P1-alpha', 'AAAA')
    second = Protein('legacy-P1-beta', 'CCCC')
    proteins = [first, second]
    assert match_protein('P1', build_id_index(proteins), proteins) is None


def test_unique_substring_fallback_remains_compatible():
    protein = Protein('legacy-P1-alpha', 'AAAA')
    proteins = [protein]
    assert match_protein('P1', build_id_index(proteins), proteins) is protein


def test_ordered_hsps_is_independent_of_blast_output_order():
    hsp_a = SimpleNamespace(query_start=20, sbjct_start=2,
                            query='CCCCCCCCCC', match='||||||||||')
    hsp_b = SimpleNamespace(query_start=10, sbjct_start=3,
                            query='AAAAAAAAAA', match='||||||||||')
    hsp_c = SimpleNamespace(query_start=10, sbjct_start=1,
                            query='BBBBBBBBBB', match='||||||||||')

    forward = SimpleNamespace(alignments=[
        SimpleNamespace(hit_id='hit-z', hsps=[hsp_a]),
        SimpleNamespace(hit_id='hit-a', hsps=[hsp_b, hsp_c]),
    ])
    reversed_ = SimpleNamespace(alignments=[
        SimpleNamespace(hit_id='hit-a', hsps=[hsp_c, hsp_b]),
        SimpleNamespace(hit_id='hit-z', hsps=[hsp_a]),
    ])

    def signature(items):
        return [(h.query_start, h.sbjct_start, h.query, h.match) for h in items]

    assert signature(ordered_hsps(forward)) == signature(ordered_hsps(reversed_))
    assert signature(ordered_hsps(forward)) == [
        (10, 1, 'BBBBBBBBBB', '||||||||||'),
        (10, 3, 'AAAAAAAAAA', '||||||||||'),
        (20, 2, 'CCCCCCCCCC', '||||||||||'),
    ]


def test_refseq_description_cannot_become_epitope_path(tmp_path):
    raw_header = ('WP_044285928.1 ZmpA/ZmpB/ZmpC family metallo-endopeptidase '
                  '[Mesomycoplasma ovipneumoniae]')
    fasta = tmp_path / 'proteome.fasta'
    fasta.write_text(f'>{raw_header}\nMAAAAAAAAA\n')

    records, _ = quality_control(str(fasta), str(tmp_path), upload=True)
    record = records[0]
    protein = Protein(str(record.description), str(record.seq))
    protein.accession = record.accession
    protein.original_id = record.original_description

    prot_name = epitope_output_identifier(protein)
    directory = tmp_path / 'epitope' / prot_name
    output_file = directory / f'mhci_epitopes_{prot_name}.csv'

    assert prot_name == 'WP_044285928.1'
    assert directory.parent == tmp_path / 'epitope'
    assert output_file.parent == directory
    assert raw_header == protein.original_id


def test_phase4_qc_policy_and_threshold():
    cleaned, reason, _, changes = qc_sequence('acUOBJZ' + 'A' * 94 + '*')
    assert cleaned == 'ACCKNLQ' + 'A' * 94
    assert reason == 'cleaned'
    assert any(change[3] == 'terminal_stop_removed' for change in changes)
    assert qc_sequence('AA*AA')[1] == 'internal_stop'
    assert qc_sequence('X' + 'A' * 19)[0].startswith('A')  # exactly 5% is kept
    assert qc_sequence('XX' + 'A' * 18)[1] == 'excessive_unknown_residues'


def test_checkpoint_is_atomic_and_signature_keyed(tmp_path):
    source = tmp_path / 'input.fasta'
    source.write_text('>p\nAAAA\n')
    signature = compute_signature([str(source)], {'gram': 'p'}, {'tree': 'abc'})
    checkpointer = Checkpointer(str(tmp_path), signature, enabled=True)
    state = {'proteins': ['p'], 'all_proteins': ['p'], 'args': 'args'}
    checkpointer.save(3, 'adhesin', state)
    assert not (tmp_path / 'checkpoint.pkl.tmp').exists()
    resumed = Checkpointer(str(tmp_path), signature, resume=True)
    assert resumed.done(3) and resumed.restored_state() == state
    source.write_text('>p\nCCCC\n')
    changed = compute_signature([str(source)], {'gram': 'p'}, {'tree': 'abc'})
    try:
        Checkpointer(str(tmp_path), changed, resume=True)
    except ValueError as error:
        assert 'does not match' in str(error)
    else:
        raise AssertionError('mismatched checkpoint was accepted')


def test_epitope_worker_autoscaling_is_bounded():
    from code.Epitope import available_cpu_count, resolve_epitope_workers
    assert available_cpu_count() >= 1
    assert resolve_epitope_workers(0, 3) == min(available_cpu_count(), 3)
    assert resolve_epitope_workers(99, 2) == 2
