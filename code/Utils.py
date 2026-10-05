#!/usr/local/bin/python
"""NERVE useful functions and classes"""

import argparse, logging, os, re, subprocess, time


class Progress:
    def __init__(self, label, total, every=30.0, min_total=1000):
        self.label = label
        self.total = total
        self.every = every
        self.enabled = total >= min_total
        self.started = time.monotonic()
        self.next_update = every

    def update(self, done):
        if not self.enabled:
            return
        elapsed = time.monotonic() - self.started
        if elapsed < self.next_update and done < self.total:
            return
        self.next_update = elapsed + self.every
        if done >= self.total:
            logging.info('%s: %s/%s in %.0fs', self.label, done, self.total, elapsed)
            return
        rate = done / elapsed if elapsed else 0
        remaining = (self.total - done) / rate if rate else 0
        logging.info('%s: %s/%s, %.0fs elapsed, ~%.0fs left',
                     self.label, done, self.total, elapsed, remaining)


def run_command(arguments, *, cwd=None):
    """Run an external tool with intact path arguments and an actionable error."""
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


# Protein identifiers.

# Characters allowed in output identifiers.
_UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')
# Characters that would break a FASTA record boundary if they reached a header.
_HEADER_CONTROL = re.compile(r'[\r\n\x00]+')
# makeblastdb -parse_seqids (Immunity.compare_proteomes) rejects the whole database when any
# accession exceeds 30 characters, and reads a '.' only as a numeric version suffix.
_VERSIONED_ACCESSION = re.compile(r'^[^.]+(\.[0-9]+)?$')
_MAX_ACCESSION = 30


def sanitize_description(description:str) -> str:
    """Return a description that cannot break FASTA record boundaries."""
    return _HEADER_CONTROL.sub('', str(description)).replace('>', '_')


def blast_safe_accession(accession:str) -> str:
    """Return an accession makeblastdb -parse_seqids can read behind sp|.

    Without this a 53-character identifier aborted the whole conservation step with
    'BLAST Database creation error: ... the accession is too long' (reproduced 2026-09-30).
    """
    # a trailing '.<version>' must be an integer; any other dot is not a version
    if not _VERSIONED_ACCESSION.match(accession):
        accession = accession.replace('.', '_')
    return accession[:_MAX_ACCESSION].rstrip('.')


def sanitize_identifier(identifier:str) -> str:
    """Return one non-empty, filesystem-safe identifier component."""
    return _UNSAFE_ID_CHARS.sub('_', str(identifier).strip()) or 'seq'


def derive_accession(header:str) -> str:
    """Extract the accession from a FASTA header."""
    stripped = header.strip()
    token = stripped.split()[0] if stripped.split() else stripped
    fields = token.split('|')
    if len(fields) >= 3 and fields[1]:
        return fields[1]
    return token


def normalize_header(header:str, seen_accessions:set) -> tuple:
    """Return a canonical FASTA header and accession."""
    stripped = header.strip()
    parts = stripped.split(None, 1)
    description = sanitize_description(parts[1]) if len(parts) > 1 else ''
    accession = blast_safe_accession(sanitize_identifier(derive_accession(stripped)))
    # Disambiguate duplicate accessions, keeping the suffixed form blast-safe too.
    if accession in seen_accessions:
        base, suffix = accession, 2
        while accession in seen_accessions:
            tail = f'_{suffix}'
            accession = blast_safe_accession(base[:_MAX_ACCESSION - len(tail)] + tail)
            suffix += 1
    seen_accessions.add(accession)
    canonical = f'sp|{accession}|{accession}'
    if description:
        canonical = f'{canonical} {description}'
    return canonical, accession


def build_id_index(list_of_proteins) -> dict:
    """Index proteins by external identifiers."""
    index = {}
    for protein in list_of_proteins:
        tokens = protein.id.split()
        index.setdefault(tokens[0] if tokens else protein.id, protein)
        if protein.accession:
            index.setdefault(protein.accession, protein)
    return index


def match_protein(external_id:str, index:dict, list_of_proteins=None):
    """Resolve an external identifier to its protein."""
    tokens = external_id.split()
    token = tokens[0] if tokens else external_id
    protein = index.get(token)
    if protein is None:
        protein = index.get(derive_accession(token))
    if protein is not None:
        return protein
    if list_of_proteins is not None:
        matches = [candidate for candidate in list_of_proteins if token in candidate.id]
        if len(matches) == 1:
            logging.warning(f'Identifier "{external_id}" resolved to "{matches[0].id}" by '
                            f'unique substring fallback; exact lookup found no match.')
            return matches[0]
        if len(matches) > 1:
            logging.error(f'Identifier "{external_id}" is ambiguous: exact lookup failed and '
                          f'{len(matches)} proteins contain "{token}"; record skipped.')
    return None


def epitope_output_identifier(protein) -> str:
    """Return the safe path component used for one protein's epitope files."""
    raw = protein.accession if protein.accession is not None else derive_accession(protein.id)
    return sanitize_identifier(raw)
