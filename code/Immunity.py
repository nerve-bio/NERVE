#!/usr/local/bin/python
"""Autoimmunity, mouse and conservation modules"""

import logging, os, tempfile
from Bio.Blast import NCBIXML 
import pandas as pd
from code import Protein
from code.Utils import build_id_index, match_protein, run_command
import shutil


def run_blast(query, database, e_value, output, threads):
    """Run BLAST only against a fresh XML destination."""
    if os.path.exists(output):
        os.remove(output)
    run_command(['blastp', '-query', query, '-db', database, '-evalue', str(e_value),
                 '-outfmt', '5', '-out', output, '-num_threads', str(threads)])
    if not os.path.isfile(output) or os.path.getsize(output) == 0:
        raise RuntimeError(f'blastp did not produce {output}')


def parse_blast_records(path):
    """Yield records from one or more concatenated BLAST XML documents."""
    document = tempfile.SpooledTemporaryFile(mode='w+t', max_size=8 * 1024 * 1024)
    document_count = 0
    try:
        with open(path) as source:
            for line in source:
                if line.lstrip().startswith('<?xml') and document.tell():
                    document.seek(0)
                    yield from NCBIXML.parse(document)
                    document.close()
                    document = tempfile.SpooledTemporaryFile(mode='w+t', max_size=8 * 1024 * 1024)
                    document_count += 1
                document.write(line)
        if document.tell():
            document.seek(0)
            yield from NCBIXML.parse(document)
            document_count += 1
        logging.debug('%s contains %s BLAST XML document(s)',
                      os.path.basename(path), document_count)
    finally:
        document.close()


def ordered_hsps(record):
    """Yield HSPs in a fixed order."""
    hsps = []
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            hsps.append((str(getattr(alignment, 'hit_id', '')),
                         int(getattr(hsp, 'query_start', 0) or 0),
                         int(getattr(hsp, 'sbjct_start', 0) or 0),
                         str(hsp.query), str(hsp.match), hsp))
    hsps.sort(key=lambda item: item[:5])
    return [item[5] for item in hsps]


def build_kmer_index(peptides, k:int) -> dict:
    """Index peptides by their k-length substrings."""
    index = {}
    for peptide in peptides:
        # Skip missing peptides.
        if not isinstance(peptide, str) or len(peptide) < k:
            continue
        for kmer in {peptide[i:i + k] for i in range(len(peptide) - k + 1)}:
            index.setdefault(kmer, []).append(peptide)
    return index


def compare_with_mhcpep(list_of_proteins, mhcpep_peptides, shared_attr:str, out_attr:str) -> None:
    """Add matching MHC peptides to each protein."""
    peptides = list(mhcpep_peptides)
    indexes = {}
    for protein in list_of_proteins:
        for seq in getattr(protein, shared_attr):
            window = seq['query'][seq['start_pos']:seq['start_pos'] + len(seq['match'])]
            if len(window) not in indexes:
                indexes[len(window)] = build_kmer_index(peptides, len(window))
            for peptide in indexes[len(window)].get(window, ()):
                getattr(protein, out_attr).extend(Protein.Protein.peptide_comparison(seq, peptide))

def autoimmunity(list_of_proteins, proteome1, working_dir, NERVE_dir, e_value, minlength, mismatch, substitution, threads=1)->list:
    """Performs research of human immunogenig peptides"""
    run_blast(proteome1, os.path.join(NERVE_dir, 'database/sapiens_database/sapiens'),
              e_value, os.path.join(working_dir, 'sapiens.xml'), threads)
    index = build_id_index(list_of_proteins)
    #logging.debug("Warning: you can find a sapiens.xml file on your working directory which is the outputs of the autoimmunity module.\nDo not delete during the computation!\nAfter the computation it will be deleted in order to avoid future collisions.")
    # for each result in the .xml file...
    #outfile = open(os.path.join(working_dir, 'autoimmunity_raw_output.txt'), 'w')
    for record in parse_blast_records(os.path.join(working_dir, "sapiens.xml")):
        query_name = record.query.split(' ')[0] # take only the query name
        # Match the query protein.
        tmp_protein = match_protein(query_name, index, list_of_proteins)
        if tmp_protein is None:
            logging.warning(f'blastp query "{query_name}" matched no input protein; record skipped.')
            continue
        # for each effective alignment between the tmp candidate and the human proteome
        for hsp in ordered_hsps(record): # collect all the interesting peptides
            tmp_protein.list_of_shared_human_peps += Protein.Protein.hsp_match_parser(hsp.match,\
                                                                                      hsp.query,\
                                                                                      parsing_window_size=minlength,\
                                                                                      max_sub=substitution,\
                                                                                      max_mismatch=mismatch)
        # print out the peptides (if there are any)
        #if len(tmp_protein.list_of_shared_human_peps) == 0:
        #    outfile.write("\nNo immunogenic peptides for " + query_name)   
        #else:
        #    outfile.write("\nList of immunogenic peptides for " + query_name + ": " +\
        #                  str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    #outfile.close()
    os.remove(os.path.join(working_dir, "sapiens.xml")) # delete after the computation
    
    # sum peptides
    logging.debug('Run sum of peptides')
    for p in list_of_proteins:
        #p.sapiens_peptides_sum=0
        score = 0
        if len(p.list_of_shared_human_peps) > 0:
            prev_match = p.list_of_shared_human_peps[0]['match']
            score = len(prev_match)
            for pept in p.list_of_shared_human_peps[1:]:
                tmp_match = pept['match']
                if tmp_match[:len(tmp_match)-1] == prev_match[1:]:
                    score += 1
                else:
                    score += len(tmp_match)
                prev_match = tmp_match
        p.sapiens_peptides_sum = score/p.length
        
    # store peptides from comparison with human recognized bacterial mhcpep
    mhcpep = pd.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_sapiens.csv"), skipinitialspace=True)
    compare_with_mhcpep(list_of_proteins, mhcpep['Epitope.2'],
                        'list_of_shared_human_peps',
                        'list_of_peptides_from_comparison_with_mhcpep_sapiens')
    return list_of_proteins
            
def mouse(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, minlength, substitution, mismatch, threads=1)->list:
    """Module to run mouse immunity check"""
    run_blast(proteome1, os.path.join(NERVE_dir, 'database/mouse_database/mouse'),
              e_value, os.path.join(working_dir, 'mouse.xml'), threads)
    index = build_id_index(list_of_proteins)
    #outfile = open(os.path.join(working_dir, 'mouse_immunity_raw_output.txt'), 'w')
    for record in parse_blast_records(os.path.join(working_dir, "mouse.xml")):
        query_name = record.query.split(' ')[0]
        # Match the query protein.
        tmp_protein = match_protein(query_name, index, list_of_proteins)
        if tmp_protein is None:
            logging.warning(f'blastp query "{query_name}" matched no input protein; record skipped.')
            continue
        for hsp in ordered_hsps(record):
            tmp_protein.list_of_shared_mouse_peps += Protein.Protein.hsp_match_parser(hsp.match, hsp.query, parsing_window_size=minlength, max_sub=substitution, max_mismatch=mismatch )
        # print out the peptides (if there are any)
        #if len(tmp_protein.list_of_shared_human_peps) == 0:
        #    outfile.write("\nNo immunogenic peptides for " + query_name)   
        #else:
        #    outfile.write("\nList of immunogenic peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    #outfile.close()
    os.remove(os.path.join(working_dir, "mouse.xml")) # delete after the computation
    
    # store peptides from comparison with mouse recognized bacterial mhcpep
    mhcpep = pd.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_mouse.csv"), skipinitialspace=True)
    compare_with_mhcpep(list_of_proteins, mhcpep['Epitope.2'],
                        'list_of_shared_mouse_peps',
                        'list_of_peptides_from_comparison_with_mhcpep_mouse')
    # sum peptides
    logging.debug('Run sum of peptides')
    for p in list_of_proteins:
        score = 0
        if len(p.list_of_shared_mouse_peps) > 0:
            prev_match = p.list_of_shared_mouse_peps[0]['match']
            score = len(prev_match)
            for pept in p.list_of_shared_mouse_peps[1:]:
                tmp_match = pept['match']
                if tmp_match[:len(tmp_match)-1] == prev_match[1:]:
                    score += 1
                else:
                    score += len(tmp_match)
                prev_match = tmp_match
        p.mouse_peptides_sum = score/p.length
    return list_of_proteins

def conservation(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, proteome2, minlength,
                 substitution, mismatch, threads=1)->list:
    """Runs conservation evaluation between proteome1 and proteome2"""
    database_dir = os.path.join(working_dir, 'compare_proteome')
    database = os.path.join(database_dir, 'compare_proteome')
    if os.path.exists(database_dir):
        shutil.rmtree(database_dir)
    os.makedirs(database_dir)
    run_command(['makeblastdb', '-in', proteome2, '-dbtype', 'prot',
                 '-parse_seqids', '-out', database])
    if not any(os.path.isfile(database + suffix) for suffix in ('.pin', '.pdb')):
        raise RuntimeError(f'makeblastdb did not produce an index at {database}')
    run_blast(proteome1, database, e_value,
              os.path.join(working_dir, 'comparison.xml'), threads)

    index = build_id_index(list_of_proteins)
    outfile=open(os.path.join(working_dir, 'conservation_raw_output.txt'), 'w')
    for record in parse_blast_records(os.path.join(working_dir, "comparison.xml")):
        query_name = record.query.split(' ')[0]

        tmp_protein = match_protein(query_name, index, list_of_proteins)
        if tmp_protein is None:
            logging.warning(f'blastp query "{query_name}" matched no input protein; record skipped.')
            continue
        for hsp in ordered_hsps(record):
            #if hsp.score > max_score: max_score = hsp.score
            tmp_protein.list_of_shared_conserv_proteome_peps += Protein.Protein.hsp_match_parser(hsp.match,
                                                                                                 hsp.query,
                                                                                                 parsing_window_size=minlength,
                                                                                                 max_sub=substitution,
                                                                                                 max_mismatch=mismatch)
        # Store conservation peptide matches.
        if len(tmp_protein.list_of_shared_conserv_proteome_peps) == 0:
            outfile.write("\nNo shared peptides for " + query_name)
        else:
            outfile.write("\nList of shared peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_conserv_proteome_peps]))
    outfile.close()                                                                                              
    # sum peptides
    logging.debug('Run sum of peptides')
    for p in list_of_proteins:
            #p.conservation_score = 0
            score = 0
            if len(p.list_of_shared_conserv_proteome_peps) > 0:
                prev_match = p.list_of_shared_conserv_proteome_peps[0]['match']
                score = len(prev_match)
                for pept in p.list_of_shared_conserv_proteome_peps[1:]:
                    tmp_match = pept['match']
                    if tmp_match[:len(tmp_match)-1] == prev_match[1:]:
                        score += 1
                    else:
                        score += len(tmp_match)
                    prev_match = tmp_match
            p.conservation_score = score/p.length

    os.remove(os.path.join(working_dir, "comparison.xml")) # delete after the computation
    if os.path.isdir(os.path.join(working_dir, "compare_proteome")):
        shutil.rmtree(os.path.join(working_dir, "compare_proteome"))
    return list_of_proteins
