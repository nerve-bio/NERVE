#!/usr/local/bin/python
"""Quality control module"""

import os, logging, re, requests
from Bio import SeqIO
from Bio.Seq import Seq
from code.Utils import dir_path, normalize_header

class protein_element:
    """Class to handle fasta file elements similarly to the biopython fasta file parser"""
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq


CANONICAL_AA = set('ACDEFGHIKLMNPQRSTVWY')
AA_SUBSTITUTIONS = {
    'U': 'C',  # selenocysteine
    'O': 'K',  # pyrrolysine
    'B': 'N',  # Asx -> neutral amide
    'Z': 'Q',  # Glx -> neutral amide
    'J': 'L',  # Xle
}
UNKNOWN_SUBSTITUTE = 'A'
DEFAULT_MAX_UNKNOWN_FRACTION = 0.05


def qc_sequence(sequence, max_unknown_fraction=DEFAULT_MAX_UNKNOWN_FRACTION):
    """Normalize one sequence and return (cleaned, reason, message, changes)."""
    if not 0 <= max_unknown_fraction <= 1:
        raise ValueError('max_unknown_fraction must be between 0 and 1')
    changes = []
    source = str(sequence).upper()
    cleaned = re.sub(r'[\s.\-\d]', '', source)
    if len(cleaned) != len(source):
        changes.append((0, f'{len(source) - len(cleaned)} char(s)', '',
                        'gaps_whitespace_digits_removed'))
    if not cleaned:
        return None, 'empty_sequence', 'sequence is empty after normalization', changes
    if '*' in cleaned:
        body = cleaned.rstrip('*')
        trailing = len(cleaned) - len(body)
        if '*' in body:
            pos = body.index('*') + 1
            return None, 'internal_stop', f'internal stop codon at position {pos}/{len(cleaned)}', changes
        if not body:
            return None, 'empty_after_stop_removal', 'sequence contained only stop symbols', changes
        changes.append((len(body) + 1, '*' * trailing, '', 'terminal_stop_removed'))
        cleaned = body
    chars = list(cleaned)
    for index, amino_acid in enumerate(chars):
        if amino_acid in AA_SUBSTITUTIONS:
            replacement = AA_SUBSTITUTIONS[amino_acid]
            chars[index] = replacement
            changes.append((index + 1, amino_acid, replacement,
                            f'ambiguity_{amino_acid}_to_{replacement}'))
    cleaned = ''.join(chars)
    unknown = [index for index, amino_acid in enumerate(cleaned)
               if amino_acid not in CANONICAL_AA]
    if unknown and len(unknown) / len(cleaned) > max_unknown_fraction:
        fraction = len(unknown) / len(cleaned)
        return None, 'excessive_unknown_residues', (
            f'{len(unknown)}/{len(cleaned)} residues ({fraction:.1%}) are unknown, '
            f'above the {max_unknown_fraction:.1%} limit'), changes
    if unknown:
        chars = list(cleaned)
        for index in unknown:
            original = chars[index]
            chars[index] = UNKNOWN_SUBSTITUTE
            changes.append((index + 1, original, UNKNOWN_SUBSTITUTE,
                            'sparse_unknown_substituted'))
        cleaned = ''.join(chars)
    return cleaned, ('cleaned' if changes else 'ok'), '', changes

def proteome_downloader(working_dir, proteome_id, filename='input_proteome.fasta', output_dir=os.getcwd(), format_ = "fasta") -> None:
    """Downloads proteome from uniprot database into output multifasta file
    param: proteome_id: uniprot unique proteome id, not case-sensitive
    param: output_dir: output directory (default: current directory)
    param: format: uniprot API required format (default:fasta)
    param: filename: output proteome filename (default: input_proteome.fasta)
    """
    # define logging file
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)
    try: 
        url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format={format_}&query=%28proteome%3A{proteome_id}%29'
        download_from_url_to_file(url, output_dir, filename, 'First download attempt failed');
    except AssertionError: 
        logging.debug(f'28proteome is not a valid building block')
        try:
            url=f'https://rest.uniprot.org/uniprotkb/stream?format={format_}&query=%28proteome%3A{proteome_id}%29'
            download_from_url_to_file(url, output_dir, filename, 'Second download attempt failed');
        except AssertionError:
            logging.debug(f'avoiding compression is not a valid building block')
            try:
                url=f'https://rest.uniprot.org/uniparc/stream?format={format_}&query=%28upid%3A{proteome_id}%29'
                download_from_url_to_file(url, output_dir, filename, 'Third download attempt failed');
            except AssertionError as e:
                logging.debug(f'28upid is not a valid building block')
                print(f'Unable to download proteome {proteome_id} due to invalid proteome ID, Uniprot API failure or wrong path. In case of uniprot API failure provide the proteome as file.')
                raise SystemExit(e)
    return None

def download_from_url_to_file(url:str, output_dir:str, filename:str, assert_error_message:str):
    downloaded = False
    for _ in range(3):
        try:
            response = requests.get(url, stream = True, timeout = 120)
            text_file = open(os.path.join(output_dir, filename), 'wb')
            for chunk in response.iter_content(chunk_size=1024):
                text_file.write(chunk)
            # raise an AssertionError if the given proteome ID is not valid
            assert text_file.tell() > 0, assert_error_message
            text_file.close()
            downloaded = True
            break
        except:
            logging.debug(f'Failed to fetch from {url}')
            continue
    assert downloaded is True, 'Download failed after retries'

def is_fasta(filename:str):
    """Function that rise an error if the format is not .fasta, the file is empty or non existing.
    param: filename: path to fasta file"""
    if not os.path.isfile(filename):
        raise ValueError("{filename} not found.")
    with open(filename, "r") as handle:
        fasta = list(SeqIO.parse(handle, "fasta"))
        # biopython silently fails if the format is not fasta returning an empty generator
        # any() returns False if the list is empty
        if any(fasta) == True:
            return fasta
        else:
            raise ValueError(f'{filename} is not in fasta format or the input fasta file is empy.')

def quality_control(path_to_fasta:str, working_dir:str, upload=False,
                    max_unknown_fraction=DEFAULT_MAX_UNKNOWN_FRACTION) -> dir_path:
    """
    Remove sequences with non-canonical aminoacid symbols. U (Se-Cys) is substituted with C (Cys). Returns
    {working_dir}/discarded_sequences_{input_fasta_file_name}.fasta containing discarded sequences and\
    {working_dir}/cleaned_{input_fasta_file_name}.fasta with cleaned sequences
    param: path_to_fasta: full path to fasta file containing the proteome with .fasta extension. Input fasta file\
    will not be overwritten
    param: working_dir: working directory, were cleaned fasta file will be saved
    param: upload: if True the function returns a list containing the filtered sequences as protein_element objects
    output: path to the cleaned fasta file named as {working_dir}/cleaned_{input_fasta_file_name}.fasta
    """
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    
    filtered_sequences, discarded_sequences = [],[]
    # accessions assigned so far, so canonical headers stay unique
    seen_accessions = set()
    # control formatting
    fasta_list = is_fasta(path_to_fasta)
    # filename needed to create the output file
    file_name = os.path.basename(path_to_fasta)
    output_file = os.path.join(working_dir, "_".join(["cleaned",file_name]))
    output_discarded_sequences=os.path.join(working_dir, "_".join(["discarded_sequences",file_name]))
    base_name = os.path.splitext(file_name)[0]
    discard_log = os.path.join(working_dir, f'qc_discard_log_{base_name}.tsv')
    modification_log = os.path.join(working_dir, f'qc_modifications_log_{base_name}.tsv')
    discard_rows, modification_rows = [], []
    for record in fasta_list:
        new_seq, reason, message, modifications = qc_sequence(
            str(record.seq), max_unknown_fraction)
        # check name: a '>' inside a header would split the record when the
        # cleaned file is parsed again, so it goes before anything else reads it
        original_description = record.description
        if ">" in original_description:
            logging.debug(f'Found non-allowed character ">" in sequence name:\n{original_description}\nSubstituting with "_"')
            original_description = original_description.replace(">", "_")
        # rewrite to a canonical, unique, filesystem-safe header so that every
        # downstream tool echoes back an identifier we can resolve exactly
        canonical_description, accession = normalize_header(original_description, seen_accessions)
        if canonical_description != original_description:
            logging.debug(f'Header normalized:\n{original_description}\n-> {canonical_description}')
        record.description = canonical_description
        # keep the user's header for reporting; the pipeline matches on the
        # canonical one but the output tables should stay readable
        record.original_description = original_description
        record.accession = accession
        if new_seq is not None:
            record.seq = Seq(new_seq)
        
        if new_seq is not None:
            filtered_sequences.append(record)
            for position, original, substituted, rule in modifications:
                modification_rows.append((accession, position, original, substituted, rule))
        else:
            discarded_sequences.append(record)
            discard_rows.append((accession, reason, message))
            logging.debug(f'DISCARD {record.description} [{reason}] {message}')
    if discard_rows:
        with open(discard_log, 'w') as handle:
            handle.write('accession\treason_code\tmessage\n')
            for row in discard_rows:
                handle.write('\t'.join(map(str, row)) + '\n')
    if modification_rows:
        with open(modification_log, 'w') as handle:
            handle.write('accession\tposition\toriginal\tsubstituted\trule\n')
            for row in modification_rows:
                handle.write('\t'.join(map(str, row)) + '\n')
    summary = (f'QC of {file_name}: {len(filtered_sequences)}/{len(fasta_list)} sequences kept, '
               f'{len(discarded_sequences)} discarded, '
               f'{len(set(row[0] for row in modification_rows))} modified.')
    print(summary)
    logging.info(summary)
    # output filtered overwriting input fasta file
    filename = open(output_file, 'w')
    for sequence in filtered_sequences:
        filename.write(f'>{str(sequence.description)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    #SeqIO.write(filtered_sequences, filename, "fasta")
    filename.close()
    # output discarded sequences
    filename = open(output_discarded_sequences, 'w')
    for sequence in discarded_sequences:
        filename.write(f'>{str(sequence.description)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    #SeqIO.write(discarded_sequences, filename, "fasta")
    filename.close()
    # repath output_file to be absolute
    output_file=os.path.join('/', os.path.relpath(output_file, start = '/'))
    if upload==True:
        # return filtered sequnces and the new path
        return filtered_sequences, output_file
    # return the path to the cleaned file
    return output_file
