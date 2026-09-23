#!/usr/local/bin/python
"""Quality control module"""

import os, logging, re, requests, tempfile
from contextlib import closing
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

def proteome_downloader(working_dir, proteome_id, filename='input_proteome.fasta',
                       output_dir=None, format_='fasta') -> None:
    """Download a UniProt proteome, publishing only a complete FASTA."""
    if format_ != 'fasta':
        raise ValueError('Only FASTA downloads are supported')
    target = filename if os.path.isabs(filename) else os.path.join(output_dir or working_dir, filename)
    urls = (
        f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28proteome%3A{proteome_id}%29',
        f'https://rest.uniprot.org/uniparc/stream?format=fasta&query=%28upid%3A{proteome_id}%29',
    )
    last_error = None
    for url in urls:
        try:
            download_from_url_to_file(url, target)
            return
        except requests.HTTPError as error:
            last_error = error
            if error.response is None or error.response.status_code not in (400, 404):
                break
        except (requests.RequestException, ValueError) as error:
            last_error = error
            break
    raise RuntimeError(f'Unable to download UniProt proteome {proteome_id}: {last_error}') from last_error


def download_from_url_to_file(url, target):
    """Download atomically; retry only transient transport or server failures."""
    directory = os.path.dirname(os.path.abspath(target))
    os.makedirs(directory, exist_ok=True)
    for attempt in range(3):
        temporary = None
        try:
            with closing(requests.get(url, stream=True, timeout=(15, 120))) as response:
                response.raise_for_status()
                with tempfile.NamedTemporaryFile(mode='wb', dir=directory, prefix='.proteome-',
                                                  delete=False) as output:
                    temporary = output.name
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            output.write(chunk)
            with open(temporary, 'rb') as handle:
                if handle.read(1) != b'>':
                    raise ValueError('Downloaded content is not FASTA')
            records = 0
            with open(temporary) as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    if not record.id or not record.seq:
                        raise ValueError('Downloaded FASTA contains an empty record')
                    records += 1
            if not records:
                raise ValueError('Downloaded FASTA contains no records')
            os.replace(temporary, target)
            return
        except requests.HTTPError as error:
            if error.response is None or error.response.status_code < 500 or attempt == 2:
                raise
        except (requests.Timeout, requests.ConnectionError,
                requests.exceptions.ChunkedEncodingError):
            if attempt == 2:
                raise
        finally:
            if temporary and os.path.exists(temporary):
                os.remove(temporary)

def is_fasta(filename:str):
    """Function that rise an error if the format is not .fasta, the file is empty or non existing.
    param: filename: path to fasta file"""
    if not os.path.isfile(filename):
        raise ValueError(f'{filename} not found.')
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
    # Keep canonical accessions unique.
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
        # Remove embedded FASTA delimiters.
        original_description = record.description
        if ">" in original_description:
            logging.debug(f'Found non-allowed character ">" in sequence name:\n{original_description}\nSubstituting with "_"')
            original_description = original_description.replace(">", "_")
        # Normalize the header for downstream tools.
        canonical_description, accession = normalize_header(original_description, seen_accessions)
        if canonical_description != original_description:
            logging.debug(f'Header normalized:\n{original_description}\n-> {canonical_description}')
        record.description = canonical_description
        # Keep the original header for reporting.
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
