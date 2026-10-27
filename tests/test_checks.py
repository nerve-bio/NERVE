#!/usr/bin/python3
"""

    Tests the input checks and the external tool guards
    run with: docker run --rm -v $(pwd):/repo -w /repo --entrypoint pytest nerve:<tag> tests/test_checks.py

"""

import os
import pytest
from code.Function import deep_fri
from code.Quality_control import quality_control, validate_downloaded_proteome
from code.Utils import MAX_SEQUENCE_LENGTH, MIN_SEQUENCE_LENGTH, length_rejection, nucleotide_fraction

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
