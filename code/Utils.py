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
            print(f'  {self.label}: {done}/{self.total} in {elapsed:.0f}s', flush=True)
            return
        rate = done / elapsed if elapsed else 0
        remaining = (self.total - done) / rate if rate else 0
        print(f'  {self.label}: {done}/{self.total}  {elapsed:.0f}s elapsed, '
              f'~{remaining:.0f}s left', flush=True)


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


# Protein identifiers.

# Characters allowed in output identifiers.
_UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')


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
    description = parts[1] if len(parts) > 1 else ''
    accession = sanitize_identifier(derive_accession(stripped))
    # Disambiguate duplicate accessions.
    if accession in seen_accessions:
        base, suffix = accession, 2
        while f'{base}_{suffix}' in seen_accessions:
            suffix += 1
        accession = f'{base}_{suffix}'
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
