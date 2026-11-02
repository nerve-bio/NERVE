#!/usr/local/bin/python
"""NERVE useful functions and classes"""

import argparse, os, re, subprocess

_SAFE_IDENTIFIER = re.compile(r'[^A-Za-z0-9._-]+')
_ACCESSION_TAGS = {'ref', 'gb', 'emb', 'dbj'}
_HEADER_CONTROL = re.compile(r'[\r\n\x00]+')
_VERSIONED_ACCESSION = re.compile(r'^[^.]+(\.[0-9]+)?$')
_NCBI_CDS_ID = re.compile(r'^(?:lcl\|)?[^|]+_prot_(.+?)_[0-9]+$')
_MAX_ACCESSION = 30
# iFeature's CTDT divides by len(sequence) - 1, CTDC by len(sequence)
MIN_SEQUENCE_LENGTH = 2
# the longest protein on record is titin, near 35,000 residues
MAX_SEQUENCE_LENGTH = 50000
# a proteome is 23-33% of these symbols, a nucleotide sequence 100%
NUCLEOTIDE_SYMBOLS = frozenset('ACGTUN')
MAX_NUCLEOTIDE_FRACTION = 0.8


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

def length_rejection(sequence, minimum=MIN_SEQUENCE_LENGTH, maximum=MAX_SEQUENCE_LENGTH):
    """Return why a sequence cannot be scored, or None"""
    if len(sequence) < minimum:
        return f'being shorter than {minimum} residues'
    if len(sequence) > maximum:
        return f'being longer than {maximum} residues'
    return None

def nucleotide_fraction(sequences):
    """Fraction of residues that are also nucleotide symbols"""
    total = matched = 0
    for sequence in sequences:
        sequence = str(sequence).upper()
        total += len(sequence)
        matched += sum(1 for residue in sequence if residue in NUCLEOTIDE_SYMBOLS)
    return matched / total if total else 0.0

CANONICAL_AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
# identity is known, so the mapping is a tool concession and costs no information
DETERMINATE_SUBSTITUTIONS = {'U': 'C', 'O': 'K'}
# Asx and Glx take the amide form, Xle takes leucine; see nerve-local/QC-THREE-WAY.md
AMBIGUOUS_SUBSTITUTIONS = {'B': 'N', 'Z': 'Q', 'J': 'L'}
# glycine is near neutral on Kyte-Doolittle, alanine is +1.8 and can fake a TM helix
UNKNOWN_SUBSTITUTE = 'G'
MAX_UNKNOWN_FRACTION = 0.05
# three times the unknown budget, from the B/Z distribution of 1,596,204 records
MAX_AMBIGUITY_FRACTION = 0.15
_STRIPPED = re.compile(r'[\s.\-\d]')
# a code-4 proteome sits at 70-76% internal stops, real pseudogene load tops out near 19%
MIN_STOP_SHARE = 0.05
MAX_STOP_SHARE = 0.50
# a share needs records to be a share; genome-scale arms clear this floor
MIN_STOP_RECORDS = 10


def sequence_case(sequence):
    """Return 'upper', 'lower' or 'mixed' for the cased residues of a sequence"""
    letters = [character for character in str(sequence) if character.isalpha()]
    if not letters:
        return 'upper'
    lower = sum(1 for character in letters if character.islower())
    if lower == 0:
        return 'upper'
    return 'lower' if lower == len(letters) else 'mixed'


def qc_sequence(sequence, max_unknown_fraction=MAX_UNKNOWN_FRACTION,
                max_ambiguity_fraction=MAX_AMBIGUITY_FRACTION):
    """Normalize one sequence and return (cleaned, reason, message, changes)"""
    changes = []
    source = str(sequence)
    case = sequence_case(source)
    if case == 'mixed':
        lower = sum(1 for character in source if character.islower())
        changes.append((0, f'{lower} of {len(source)} lowercase', '', 'mixed_case_uppercased'))
    source = source.upper()
    cleaned = _STRIPPED.sub('', source)
    if len(cleaned) != len(source):
        changes.append((0, f'{len(source) - len(cleaned)} char(s)', '',
                        'gaps_whitespace_digits_removed'))
    if not cleaned:
        return None, 'empty_sequence', 'sequence is empty after normalization', changes
    if '*' in cleaned:
        body = cleaned.rstrip('*')
        trailing = len(cleaned) - len(body)
        if '*' in body:
            position = body.index('*') + 1
            return None, 'internal_stop', (
                f'internal stop codon at position {position}/{len(cleaned)}'), changes
        if not body:
            return None, 'empty_after_stop_removal', 'sequence contained only stop symbols', changes
        changes.append((len(body) + 1, '*' * trailing, '', 'terminal_stop_removed'))
        cleaned = body
    # length is measured here so a terminator cannot carry a record over the floor
    rejection = length_rejection(cleaned)
    if rejection is not None:
        return None, 'unscorable_length', rejection, changes
    # '?' claims a position of unknown content, so it joins the unknown budget
    if '?' in cleaned:
        changes.append((cleaned.index('?') + 1, '?', 'X', 'missing_marked_unknown'))
        cleaned = cleaned.replace('?', 'X')
    # anything left outside the alphabet is not a residue code and is not filled
    symbol = next((character for character in cleaned if not character.isalpha()), None)
    if symbol is not None:
        return None, 'non_residue_symbol', (
            f'"{symbol}" at position {cleaned.index(symbol) + 1} is not a residue symbol'), changes
    chars = list(cleaned)
    for index, amino_acid in enumerate(chars):
        if amino_acid in DETERMINATE_SUBSTITUTIONS:
            chars[index] = DETERMINATE_SUBSTITUTIONS[amino_acid]
            changes.append((index + 1, amino_acid, chars[index],
                            f'determinate_{amino_acid}_to_{chars[index]}'))
    ambiguous = [index for index, amino_acid in enumerate(chars)
                 if amino_acid in AMBIGUOUS_SUBSTITUTIONS]
    if ambiguous and len(ambiguous) / len(chars) > max_ambiguity_fraction:
        return None, 'excessive_ambiguity_residues', (
            f'{len(ambiguous)}/{len(chars)} residues ({len(ambiguous) / len(chars):.1%}) are '
            f'ambiguous, above the {max_ambiguity_fraction:.1%} limit'), changes
    for index in ambiguous:
        amino_acid = chars[index]
        chars[index] = AMBIGUOUS_SUBSTITUTIONS[amino_acid]
        changes.append((index + 1, amino_acid, chars[index],
                        f'ambiguity_{amino_acid}_to_{chars[index]}'))
    unknown = [index for index, amino_acid in enumerate(chars) if amino_acid not in CANONICAL_AA]
    if unknown and len(unknown) / len(chars) > max_unknown_fraction:
        return None, 'excessive_unknown_residues', (
            f'{len(unknown)}/{len(chars)} residues ({len(unknown) / len(chars):.1%}) are unknown, '
            f'above the {max_unknown_fraction:.1%} limit'), changes
    for index in unknown:
        changes.append((index + 1, chars[index], UNKNOWN_SUBSTITUTE, 'sparse_unknown_substituted'))
        chars[index] = UNKNOWN_SUBSTITUTE
    cleaned = ''.join(chars)
    return cleaned, ('cleaned' if changes else 'ok'), '', changes

def run_command(arguments, *, cwd=None):
    """Run an external tool, raising with its exit code and stderr"""
    result = subprocess.run(arguments, cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if result.returncode:
        detail = result.stderr.strip()[-1000:]
        raise RuntimeError(f'{arguments[0]} exited {result.returncode}: {detail}')
    return result.stdout, result.stderr

def dir_path(path:str) -> str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path
