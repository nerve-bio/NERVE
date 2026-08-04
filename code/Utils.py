#!/usr/local/bin/python
"""NERVE useful functions and classes"""

import argparse, logging, os, re, subprocess

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


# --- protein identity -------------------------------------------------------
# Every external tool in the pipeline (iFeature, psortb, blastp, DeepFRI) echoes
# back an identifier derived from the FASTA header, and each module has to map
# that string onto the Protein it came from. Doing that with `if tool_id in
# protein.id` is wrong twice over: it is a scan of the whole proteome per lookup,
# and it matches whenever one identifier is a substring of another, so a short id
# silently binds to the first longer protein that contains it.
#
# Instead headers are rewritten to a canonical, unique, filesystem-safe form
# during quality control (normalize_header), and resolved here by exact lookup.
# Uniqueness is what makes the exact lookup provably unambiguous rather than
# merely better than substring matching.

# anything outside this set is replaced: the accession is also used as a
# directory name for the per-protein epitope output
_UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')


def derive_accession(header:str) -> str:
    """Extract the accession from a FASTA header (without the leading '>').

    'db|ACC|NAME' style headers (UniProt, and the NCBI variants that carry more
    than three fields) yield the second '|'-delimited field; anything else -- a
    bare or in-house identifier -- yields the first whitespace-delimited token.
    param: header: FASTA header without the leading '>'
    output: raw, un-sanitized accession string
    """
    stripped = header.strip()
    token = stripped.split()[0] if stripped.split() else stripped
    fields = token.split('|')
    if len(fields) >= 3 and fields[1]:
        return fields[1]
    return token


def normalize_header(header:str, seen_accessions:set) -> tuple:
    """Rewrite a FASTA header to canonical 'sp|ACC|ACC description' form.

    ACC is sanitized to be safe as a filename and made unique across the
    proteome, so that identifier matching and the per-protein epitope output
    directories both work whatever convention the input database follows.
    param: header: FASTA header without the leading '>'
    param: seen_accessions: accessions already assigned; mutated in place
    output: (canonical_header, accession)
    """
    stripped = header.strip()
    parts = stripped.split(None, 1)
    description = parts[1] if len(parts) > 1 else ''
    accession = _UNSAFE_ID_CHARS.sub('_', derive_accession(stripped)) or 'seq'
    # disambiguate collisions so the accession identifies exactly one protein
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
    """Index proteins by the identifiers external tools can emit for them.

    Two aliases are registered per protein: the leading token of its id (the
    canonical 'sp|ACC|ACC') and the bare accession, since tools differ in how
    much of the header they echo back. normalize_header makes both unique, so
    the mapping is unambiguous; the first protein registered wins if it is not,
    which keeps the result deterministic on un-normalized input.
    param: list_of_proteins: list of Protein objects
    output: dict mapping identifier alias -> Protein
    """
    index = {}
    for protein in list_of_proteins:
        tokens = protein.id.split()
        index.setdefault(tokens[0] if tokens else protein.id, protein)
        if protein.accession:
            index.setdefault(protein.accession, protein)
    return index


def match_protein(external_id:str, index:dict, list_of_proteins=None):
    """Resolve an identifier emitted by an external tool to its Protein.

    Tries the leading token, then the accession derived from it. Both are exact
    dictionary lookups. A miss cannot happen for normalized headers, so the
    substring fallback exists only to keep un-normalized input working and warns
    when it fires.
    param: external_id: identifier string as emitted by an external tool
    param: index: mapping produced by build_id_index
    param: list_of_proteins: enables the substring fallback when supplied
    output: the matching Protein, or None
    """
    tokens = external_id.split()
    token = tokens[0] if tokens else external_id
    protein = index.get(token)
    if protein is None:
        protein = index.get(derive_accession(token))
    if protein is not None:
        return protein
    if list_of_proteins is not None:
        for candidate in list_of_proteins:
            if token in candidate.id:
                logging.warning(f'Identifier "{external_id}" resolved to "{candidate.id}" by '
                                f'substring fallback; exact lookup found no match.')
                return candidate
    return None
