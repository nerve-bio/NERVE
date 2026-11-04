#!/usr/bin/python3
"""

    Tests the input checks and the external tool guards
    run with: docker run --rm -v $(pwd):/repo -w /repo --entrypoint pytest nerve:<tag> tests/test_checks.py

"""

import os
import pytest
import shutil
from code.Function import deep_fri
from code.Quality_control import quality_control, validate_downloaded_proteome
from code.Utils import MAX_AMBIGUITY_FRACTION, MAX_SEQUENCE_LENGTH, MIN_SEQUENCE_LENGTH, UNKNOWN_SUBSTITUTE, distinct_input, length_rejection, nucleotide_fraction, qc_sequence, sequence_case

GOOD = '>sp|A1|A1 one\nMKTAYIAKQR\n>sp|A2|A2 two\nMQIFVKTLTG\n'
FAILING_PREDICT = 'import sys\nsys.stderr.write("CUDA out of memory\\n")\nsys.exit(1)\n'

def write(directory, name, text):
    """Writes text to directory/name and returns the path"""
    path = os.path.join(str(directory), name)
    open(path, 'w').write(text)
    return path

def failing_deepfri(tmp_path):
    """Builds a DeepFri directory whose predict.py exits non zero"""
    deepfri_dir = os.path.join(str(tmp_path), 'DeepFri')
    os.mkdir(deepfri_dir)
    write(deepfri_dir, 'predict.py', FAILING_PREDICT)
    return deepfri_dir

def test_length_rejection_bounds():
    """Both bounds are inclusive"""
    assert length_rejection('M' * MIN_SEQUENCE_LENGTH) is None
    assert length_rejection('M' * MAX_SEQUENCE_LENGTH) is None
    assert 'shorter' in length_rejection('M' * (MIN_SEQUENCE_LENGTH - 1))
    assert 'longer' in length_rejection('M' * (MAX_SEQUENCE_LENGTH + 1))

def test_nucleotide_fraction():
    """Counts ACGTUN in either case over every sequence"""
    assert nucleotide_fraction(['ACGT']) == 1.0
    assert nucleotide_fraction(['acgt']) == 1.0
    assert nucleotide_fraction(['MKLV']) == 0.0
    assert nucleotide_fraction(['ACGT', 'MKLV']) == 0.5
    assert nucleotide_fraction([]) == 0.0

def test_quality_control_discards_unscorable_lengths(tmp_path):
    """An empty, a one residue and an over cap record leave only the two good ones"""
    fasta = write(tmp_path, 'in.fasta', GOOD + '>sp|EMPTY|EMPTY\n\n>sp|SHORT|SHORT\nM\n'
                  + f'>sp|LONG|LONG\n{"M" * (MAX_SEQUENCE_LENGTH + 1)}\n')
    filtered_sequences, cleaned = quality_control(fasta, str(tmp_path), upload=True)
    assert [sequence.accession for sequence in filtered_sequences] == ['A1', 'A2']
    assert open(cleaned).read().count('>') == 2
    assert open(os.path.join(str(tmp_path), 'discarded_sequences_in.fasta')).read().count('>') == 3

def test_quality_control_reports_nucleotide_input(tmp_path, capsys):
    """A nucleotide file is reported, a protein file is not"""
    quality_control(write(tmp_path, 'dna.fasta', '>sp|G1|G1\n' + 'ACGT' * 100 + '\n'), str(tmp_path))
    assert 'may be a nucleotide sequence' in capsys.readouterr().out
    quality_control(write(tmp_path, 'prot.fasta', GOOD), str(tmp_path))
    assert 'may be a nucleotide sequence' not in capsys.readouterr().out

def test_validate_downloaded_proteome_accepts_a_whole_response(tmp_path):
    """A complete response passes"""
    validate_downloaded_proteome(write(tmp_path, 'in.fasta', GOOD), 'UP000001584')

def test_validate_downloaded_proteome_rejects_a_cut_response(tmp_path):
    """A response cut after a header line leaves that record empty"""
    path = write(tmp_path, 'in.fasta', GOOD + '>sp|A3|A3 three\n')
    with pytest.raises(ValueError, match='empty record'):
        validate_downloaded_proteome(path, 'UP000001584')

def test_validate_downloaded_proteome_rejects_an_error_page(tmp_path):
    """An API error page is not fasta at all"""
    path = write(tmp_path, 'in.fasta', '<html><body>503 Service Unavailable</body></html>\n')
    with pytest.raises(ValueError, match='not in fasta format'):
        validate_downloaded_proteome(path, 'UP000001584')

def test_deep_fri_reports_the_exit_code(tmp_path):
    """A failed predict.py is raised with its code, not as a missing output file"""
    with pytest.raises(RuntimeError, match='exited 1'):
        deep_fri(write(tmp_path, 'in.fasta', GOOD), failing_deepfri(tmp_path), str(tmp_path) + '/')

def test_deep_fri_does_not_parse_a_stale_file(tmp_path):
    """A previous run's predictions are removed before predict.py runs"""
    write(tmp_path, '_MF_predictions.csv',
          'Protein,GO_term/EC_number,Score,GO_term/EC_number name\nsp|STALE|STALE,GO:1,0.9,stale\n')
    with pytest.raises(RuntimeError):
        deep_fri(write(tmp_path, 'in.fasta', GOOD), failing_deepfri(tmp_path), str(tmp_path) + '/')
    assert not os.path.exists(os.path.join(str(tmp_path), '_MF_predictions.csv'))

def test_determinate_substitutions_are_exempt():
    """U and O are resolved and never counted against a budget"""
    cleaned, reason, _, _ = qc_sequence('U' * 50)
    assert cleaned == 'C' * 50 and reason == 'cleaned'
    cleaned, reason, _, _ = qc_sequence('O' * 50)
    assert cleaned == 'K' * 50 and reason == 'cleaned'

def test_ambiguity_budget_bounds_b_z_and_j():
    """B, Z and J share their own budget, which the donor left unbounded"""
    assert qc_sequence('B' * 50)[1] == 'excessive_ambiguity_residues'
    assert qc_sequence('J' * 50)[1] == 'excessive_ambiguity_residues'
    under = 'B' * 7 + 'A' * 43
    assert len(under) * MAX_AMBIGUITY_FRACTION >= 7
    assert qc_sequence(under)[0] == 'N' * 7 + 'A' * 43

def test_unknown_budget_counts_x_alone():
    """X keeps the donor's 5% limit and is filled, not resolved"""
    assert qc_sequence('X' + 'A' * 49)[0] == UNKNOWN_SUBSTITUTE + 'A' * 49
    assert qc_sequence('XXX' + 'A' * 47)[1] == 'excessive_unknown_residues'

def test_malformed_symbols_are_refused_but_question_marks_are_not():
    """A pipe is not a residue, a question mark is a position of unknown content"""
    assert qc_sequence('|' + 'A' * 99)[1] == 'non_residue_symbol'
    assert qc_sequence('?' + 'A' * 99)[0] == UNKNOWN_SUBSTITUTE + 'A' * 99

def test_stop_codons_strip_at_the_end_and_discard_inside():
    """A terminator is legitimate, an internal stop means the record does not exist"""
    assert qc_sequence('MKTAYIAKQR***')[0] == 'MKTAYIAKQR'
    assert qc_sequence('MKTA*YIAKQR')[1] == 'internal_stop'

def test_length_is_measured_after_the_stop_strip():
    """M* is a one residue protein, not a two residue one"""
    assert qc_sequence('M*')[1] == 'unscorable_length'
    assert qc_sequence('MA')[1] == 'ok'

def test_case_is_normalized_and_mixed_case_is_disclosed():
    """Uppercasing rescues a lowercase file; only a mask is worth reporting"""
    assert qc_sequence('mktayiakqr')[0] == 'MKTAYIAKQR'
    assert sequence_case('MKTA') == 'upper'
    assert sequence_case('mkta') == 'lower'
    assert sequence_case('MKta') == 'mixed'
    rules = [change[3] for change in qc_sequence('MKta')[3]]
    assert 'mixed_case_uppercased' in rules

def test_quality_control_reports_the_internal_stop_share(tmp_path, capsys):
    """Below the floor nothing is said, above it the translation table is named"""
    def proteome(carrying):
        return write(tmp_path, 'in.fasta', ''.join(
            '>sp|S%d|S%d\n%s\n' % (index, index, 'MKTA*YIAKQR' if index < carrying else 'MKTAYIAKQR')
            for index in range(100)))
    quality_control(proteome(3), str(tmp_path))
    assert 'internal stop' not in capsys.readouterr().out
    quality_control(proteome(20), str(tmp_path))
    reported = capsys.readouterr().out
    assert '20 of 100' in reported and 'translation table' not in reported
    quality_control(proteome(75), str(tmp_path))
    assert 'translation table' in capsys.readouterr().out

def test_quality_control_writes_the_audit_tables(tmp_path):
    """Every discard and every modification is recorded per record"""
    fasta = write(tmp_path, 'in.fasta', GOOD + '>sp|SEL|SEL\nMKTUYIAKQR\n>sp|BAD|BAD\nMKTA*YIAKQR\n')
    quality_control(fasta, str(tmp_path))
    discards = open(os.path.join(str(tmp_path), 'qc_discard_log_in.tsv')).read()
    assert 'BAD' in discards and 'internal_stop' in discards
    modifications = open(os.path.join(str(tmp_path), 'qc_modifications_log_in.tsv')).read()
    assert 'SEL' in modifications and 'determinate_U_to_C' in modifications

def test_quality_control_renames_repeated_identifiers(tmp_path, capsys):
    """The pipeline renames rather than refuses, and says how many"""
    fasta = write(tmp_path, 'in.fasta', '>sp|A|A one\nMKTAYIAKQR\n>sp|A|A two\nMQIFVKTLTG\n')
    filtered_sequences, _ = quality_control(fasta, str(tmp_path), upload=True)
    assert [sequence.accession for sequence in filtered_sequences] == ['A', 'A_2']
    assert '1 repeated identifier' in capsys.readouterr().out

def test_distinct_input_leaves_unambiguous_paths_alone(tmp_path):
    """Different filenames write different outputs, so nothing is staged"""
    first = write(tmp_path, 'strainA.fasta', GOOD)
    second = write(tmp_path, 'strainB.fasta', GOOD)
    assert distinct_input(second, first, str(tmp_path), 'proteome2') == second
    assert distinct_input(second, None, str(tmp_path), 'proteome2') == second
    assert distinct_input(first, first, str(tmp_path), 'proteome2') == first

def test_distinct_input_stages_a_shared_filename(tmp_path):
    """Two proteomes named alike in different directories would overwrite each other"""
    one, two = os.path.join(str(tmp_path), 'a'), os.path.join(str(tmp_path), 'b')
    os.mkdir(one)
    os.mkdir(two)
    first = write(one, 'proteome.fasta', GOOD)
    second = write(two, 'proteome.fasta', GOOD)
    staged = distinct_input(second, first, str(tmp_path), 'proteome2')
    assert staged != second
    assert os.path.basename(staged) != os.path.basename(first)

def test_distinct_input_escapes_a_colliding_staging_name(tmp_path):
    """The staging name itself must not collide with proteome1"""
    one, two = os.path.join(str(tmp_path), 'a'), os.path.join(str(tmp_path), 'b')
    os.mkdir(one)
    os.mkdir(two)
    first = write(one, 'proteome2.fasta', GOOD)
    second = write(two, 'proteome2.fasta', GOOD)
    staged = distinct_input(second, first, str(tmp_path), 'proteome2')
    assert os.path.basename(staged) != os.path.basename(first)
    assert not os.path.exists(staged)

def test_quality_control_outputs_stay_separate_after_staging(tmp_path):
    """The staged path produces its own cleaned file, so conservation gets two proteomes"""
    one, two = os.path.join(str(tmp_path), 'a'), os.path.join(str(tmp_path), 'b')
    working = os.path.join(str(tmp_path), 'wd')
    for directory in (one, two, working):
        os.mkdir(directory)
    first = write(one, 'proteome.fasta', '>sp|P1|P1 first\nMKTAYIAKQR\n')
    second = write(two, 'proteome.fasta', '>sp|P2|P2 second\nMQIFVKTLTG\n')
    staged = distinct_input(second, first, working, 'proteome2')
    shutil.copyfile(second, staged)
    _, cleaned_first = quality_control(first, working, upload=True)
    cleaned_second = quality_control(staged, working)
    assert os.path.abspath(cleaned_first) != os.path.abspath(cleaned_second)
    assert 'P1' in open(cleaned_first).read()
    assert 'P2' in open(cleaned_second).read()

def test_quality_control_accepts_a_bare_header(tmp_path):
    """A header that is only '>' splits to no fields, so the duplicate count must not index it"""
    working = os.path.join(str(tmp_path), 'wd')
    os.mkdir(working)
    source = write(str(tmp_path), 'bare.fasta', '>\nMKTAYIAKQR\n>sp|P1|P1 named\nMQIFVKTLTG\n')
    kept, cleaned = quality_control(source, working, upload=True)
    assert len(kept) == 2
    assert 'MKTAYIAKQR' in open(cleaned).read()
