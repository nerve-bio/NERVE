#!/usr/local/bin/python
"""Quality control module"""

import os, logging, requests
from Bio import SeqIO
from Bio.Seq import Seq
from code.Utils import dir_path, normalize_header, nucleotide_fraction, \
                       qc_sequence, MAX_NUCLEOTIDE_FRACTION, MAX_STOP_SHARE, \
                       MIN_STOP_RECORDS, MIN_STOP_SHARE

class protein_element:
    """Class to handle fasta file elements similarly to the biopython fasta file parser"""
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

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
    validate_downloaded_proteome(os.path.join(output_dir, filename), proteome_id)
    return None

def validate_downloaded_proteome(path:str, proteome_id:str) -> None:
    """A response cut after a header line is still valid fasta, with that record empty"""
    advice = 'Retry the job, or upload the proteome as a file.'
    try:
        fasta_list = is_fasta(path)
    except ValueError:
        raise ValueError(f'The download for {proteome_id} is not in fasta format; it may be a Uniprot error page. {advice}')
    for record in fasta_list:
        if len(record.seq) == 0:
            raise ValueError(f'Downloaded FASTA for {proteome_id} contains an empty record; the download may be truncated. {advice}')

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
def warn(message:str) -> None:
    """Report a proteome-level observation without refusing the job"""
    print(f'WARNING: {message}')
    logging.warning(message)

def write_audit(path:str, header:str, rows:list) -> None:
    """Write one per-record audit table, or nothing when there is nothing to record"""
    if not rows:
        return
    with open(path, 'w') as handle:
        handle.write(f'{header}\n')
        for row in rows:
            handle.write('\t'.join(map(str, row)) + '\n')

def stop_warning(file_name:str, discard_rows:list, records:int) -> None:
    """Internal stops are a file-level symptom, so report the share and never refuse"""
    stops = sum(1 for row in discard_rows if row[1] == 'internal_stop')
    if records < MIN_STOP_RECORDS or stops / records < MIN_STOP_SHARE:
        return
    share = stops / records
    message = f'{file_name}: {stops} of {records} records ({share:.0%}) carry an internal stop codon.'
    if share > MAX_STOP_SHARE:
        message += (' The translation table is probably wrong; genetic code 4 reads TGA as'
                    ' tryptophan. Check the table and the reading frame.')
    warn(message)

def quality_control(path_to_fasta:str, working_dir:str, upload=False) -> dir_path:
    """
    Normalize sequences and discard what cannot be scored. Returns
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
    seen_accessions, seen_tokens = set(), set()
    discard_rows, modification_rows = [], []
    duplicates = 0
    # control formatting
    fasta_list = is_fasta(path_to_fasta)
    # filename needed to create the output file
    file_name = os.path.basename(path_to_fasta)
    output_file = os.path.join(working_dir, "_".join(["cleaned",file_name]))
    output_discarded_sequences=os.path.join(working_dir, "_".join(["discarded_sequences",file_name]))
    base_name = os.path.splitext(file_name)[0]
    discard_log = os.path.join(working_dir, f'qc_discard_log_{base_name}.tsv')
    modification_log = os.path.join(working_dir, f'qc_modifications_log_{base_name}.tsv')
    # measured over the raw records, before any substitution puts C, N or G in them
    fraction = nucleotide_fraction(record.seq for record in fasta_list)
    if fraction > MAX_NUCLEOTIDE_FRACTION:
        warn(f'{file_name} is {fraction:.0%} nucleotide symbols and may be a nucleotide sequence, not a proteome.')
    for record in fasta_list:
        new_seq, reason, message, modifications = qc_sequence(str(record.seq))
        # standardize identifier and retain original name
        record.original_description = record.description
        # a header can be a bare '>', which splits to nothing; normalize_header renames those too
        fields = record.description.split(None, 1)
        token = fields[0] if fields else ''
        if token in seen_tokens:
            duplicates += 1
        seen_tokens.add(token)
        record.description, record.accession = normalize_header(record.description, seen_accessions)
        if new_seq is not None:
            record.seq = Seq(new_seq)
            filtered_sequences.append(record)
            for position, original, substituted, rule in modifications:
                modification_rows.append((record.accession, position, original, substituted, rule))
        else:
            discarded_sequences.append(record)
            discard_rows.append((record.accession, reason, message))
            logging.debug(f'Sequence {record.description} has been discarded for {reason}: {message}.')
    stop_warning(file_name, discard_rows, len(fasta_list))
    if duplicates:
        warn(f'{file_name} contains {duplicates} repeated identifier(s), renamed to keep them distinct.')
    write_audit(discard_log, 'accession\treason_code\tmessage', discard_rows)
    write_audit(modification_log, 'accession\tposition\toriginal\tsubstituted\trule', modification_rows)
    logging.info(f'QC of {file_name}: {len(filtered_sequences)}/{len(fasta_list)} sequences kept, '
                 f'{len(discarded_sequences)} discarded, '
                 f'{len(set(row[0] for row in modification_rows))} modified.')
    # output filtered overwriting input fasta file
    filename = open(output_file, 'w')
    for sequence in filtered_sequences:
        filename.write(f'>{str(sequence.description)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    filename.close()
    # output discarded sequences
    filename = open(output_discarded_sequences, 'w')
    for sequence in discarded_sequences:
        filename.write(f'>{str(sequence.description)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    filename.close()
    # repath output_file to be absolute
    output_file=os.path.join('/', os.path.relpath(output_file, start = '/'))
    if upload==True:
        # return filtered sequnces and the new path
        return filtered_sequences, output_file
    # return the path to the cleaned file
    return output_file
