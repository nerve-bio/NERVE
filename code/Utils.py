#!/usr/local/bin/python
"""NERVE useful functions and classes"""

import argparse, os, re, subprocess

_SAFE_IDENTIFIER = re.compile(r'[^A-Za-z0-9._-]+')
_ACCESSION_TAGS = {'ref', 'gb', 'emb', 'dbj'}
_HEADER_CONTROL = re.compile(r'[\r\n\x00]+')
_VERSIONED_ACCESSION = re.compile(r'^[^.]+(\.[0-9]+)?$')
_NCBI_CDS_ID = re.compile(r'^(?:lcl\|)?[^|]+_prot_(.+?)_[0-9]+$')
_MAX_ACCESSION = 30


def blast_safe_accession(accession):
    """Return an accession makeblastdb -parse_seqids can read behind sp|."""
    # a trailing .<version> must be an integer
    if not _VERSIONED_ACCESSION.match(accession):
        accession = accession.replace('.', '_')
    # 30 characters max
    return accession[:_MAX_ACCESSION].rstrip('.')


def sanitize_description(description):
    """Return a description that cannot break FASTA record boundaries."""
    return _HEADER_CONTROL.sub('', str(description)).replace('>', '_')


def sanitize_identifier(identifier):
    """Return a filesystem-safe, non-empty protein identifier."""
    safe = _SAFE_IDENTIFIER.sub('_', str(identifier).strip())
    safe = safe.strip('._-')
    return safe or 'seq'


def derive_accession(header):
    """Extract a useful accession from common FASTA identifier conventions."""
    token = str(header).lstrip('>').strip().split(None, 1)[0] if str(header).lstrip('>').strip() else ''
    parts = token.split('|')

    # NCBI CDS translations: <nucleotide>_prot_<protein>_<n>
    cds = _NCBI_CDS_ID.match(token)
    if cds:
        return sanitize_identifier(cds.group(1))

    if parts and parts[0].lower() == 'gi':
        for index, part in enumerate(parts[:-1]):
            if part.lower() in _ACCESSION_TAGS and parts[index + 1]:
                return sanitize_identifier(parts[index + 1])
    if len(parts) >= 3 and parts[0].lower() == 'gnl':
        return sanitize_identifier(parts[2])
    if len(parts) >= 2 and parts[0].lower() in {'sp', 'tr', 'ref', 'gb', 'emb', 'dbj', 'lcl', 'pdb'}:
        return sanitize_identifier(parts[1])
    if len(parts) >= 3 and parts[0].lower() in {'pir', 'prf'}:
        return sanitize_identifier(next((part for part in parts[1:] if part), token))
    if len(parts) >= 2 and parts[1]:
        return sanitize_identifier(parts[1])
    return sanitize_identifier(token)


def normalize_header(header, seen_accessions=None):
    """Canonicalize a FASTA header while retaining its human-readable description."""
    original = str(header).lstrip('>').strip()
    header_parts = original.split(None, 1)
    first = header_parts[0] if header_parts else ''
    description = sanitize_description(header_parts[1]) if len(header_parts) == 2 else ''
    accession = blast_safe_accession(derive_accession(first))

    if seen_accessions is not None:
        base = accession
        suffix = 2
        while accession in seen_accessions:
            tail = f'_{suffix}'
            accession = blast_safe_accession(base[:_MAX_ACCESSION - len(tail)] + tail)
            suffix += 1
        seen_accessions.add(accession)

    canonical = f'sp|{accession}|{accession}'
    if description:
        canonical += f' {description}'
    return canonical, accession


def _protein_aliases(protein):
    aliases = set()
    for value in (getattr(protein, 'accession', None), getattr(protein, 'id', None),
                  getattr(protein, 'original_id', None)):
        if value:
            value = str(value).lstrip('>').strip()
            aliases.add(value)
            aliases.add(value.split(None, 1)[0])
            aliases.add(derive_accession(value))
    return aliases


def build_id_index(proteins):
    """Build an exact alias index; ambiguous aliases deliberately map to None."""
    index = {}
    for protein in proteins:
        for alias in _protein_aliases(protein):
            if alias in index and index[alias] is not protein:
                index[alias] = None
            else:
                index[alias] = protein
    return index


def match_protein(identifier, index, proteins=None):
    """Match tool output exactly, with a unique-only fallback for legacy output."""
    value = str(identifier).lstrip('>').strip()
    candidates = (value, value.split(None, 1)[0], derive_accession(value))
    for candidate in candidates:
        if candidate in index and index[candidate] is not None:
            return index[candidate]

    if proteins is None:
        return None
    matches = [protein for protein in proteins
               if any(value in alias or alias in value for alias in _protein_aliases(protein))]
    unique = []
    for protein in matches:
        if all(protein is not existing for existing in unique):
            unique.append(protein)
    return unique[0] if len(unique) == 1 else None


def epitope_output_identifier(protein):
    """Return the safe accession used for epitope directories and filenames."""
    return sanitize_identifier(getattr(protein, 'accession', None) or derive_accession(protein.id))

def bashCmdMethod(bashCmd):
    """Run bash commands
    param: bashCmd: bash command to be run"""
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error
    
def dir_path(path:str) -> str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path
