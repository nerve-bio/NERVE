#!/usr/local/bin/python
"""Autoimmunity, mouse and conservation modules"""

import logging, os, subprocess
from Bio.Blast.Applications import NcbiblastpCommandline    
from Bio.Blast import NCBIXML 
import pandas as pd
from code import Protein
import shutil

def autoimmunity(list_of_proteins, proteome1, working_dir, NERVE_dir, e_value, minlength, mismatch, substitution)->list:
    """Performs research of human immunogenig peptides"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    
    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir, "database/sapiens_database/sapiens"), \
                                         evalue=e_value, outfmt=5, out=os.path.join(working_dir,"sapiens.xml")) # 5 is for xml 
    stdout, stderr = blastx_cline()
    #logging.debug("Warning: you can find a sapiens.xml file on your working directory which is the outputs of the autoimmunity module.\nDo not delete during the computation!\nAfter the computation it will be deleted in order to avoid future collisions.")
    # for each result in the .xml file...
    #outfile = open(os.path.join(working_dir, 'autoimmunity_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir,"sapiens.xml"))):
        query_name = record.query.split(' ')[0] # take only the query name 
        # take the right candidate to update
        for p in list_of_proteins:
            if query_name in p.id: # do not use query_name == p.id
                tmp_protein = p
        # for each effective alignment between the tmp candidate and the human proteome
        for alignment in record.alignments:
            for hsp in alignment.hsps: # collect all the interesting peptides
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
    #number_of_proteins = len(list_of_proteins)
    for p in list_of_proteins:
        for seq in p.list_of_shared_human_peps:
            #print(p.id)
            for pep in mhcpep['Epitope.2']:
                tmp_matches = Protein.Protein.peptide_comparison(seq, pep)
                #if len(str(tmp_matches))>3:
                #    print(tmp_matches)
                p.list_of_peptides_from_comparison_with_mhcpep_sapiens += tmp_matches
    return list_of_proteins
            
def mouse(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, minlength, substitution, mismatch)->list:
    """Module to run mouse immunity check"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir,"database/mouse_database/mouse"), evalue=e_value, outfmt=5, out=os.path.join(working_dir, "mouse.xml"))
    stdout, stderr = blastx_cline()
    #outfile = open(os.path.join(working_dir, 'mouse_immunity_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir, "mouse.xml"))):
        query_name = record.query.split(' ')[0]
        tmp_protein = list_of_proteins[0]
        # take the right protein
        for p in list_of_proteins:
            if query_name in p.id: # do not use query_name == p.id
                tmp_protein = p
        for alignment in record.alignments:
            for hsp in alignment.hsps:
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
    #number_of_proteins = len(list_of_proteins)
    for p in list_of_proteins:
        for seq in p.list_of_shared_mouse_peps:
            for pep in mhcpep['Epitope.2']:
                tmp_matches = Protein.Protein.peptide_comparison(seq, pep)
                p.list_of_peptides_from_comparison_with_mhcpep_mouse += tmp_matches
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
                 substitution, mismatch)->list:
    """Runs conservation evaluation between proteome1 and proteome2"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    
    bashCmd = f"makeblastdb -in {proteome2} -dbtype prot -parse_seqids -out {os.path.join(working_dir, 'compare_proteome/compare_proteome')}"
    process = subprocess.Popen(bashCmd.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()
    
    blastx_cline = NcbiblastpCommandline(query = proteome1, db = os.path.join(working_dir, 'compare_proteome/compare_proteome'), evalue = e_value, outfmt = 5, out=os.path.join(working_dir,"comparison.xml")) # 5 is for xml 
    stdout, stderr = blastx_cline()
    
    outfile=open(os.path.join(working_dir, 'conservation_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir,"comparison.xml"))):
        query_name = record.query.split(' ')[0] 
    
        for p in list_of_proteins:
            if query_name in p.id: # do not use p.id == query_name
                tmp_protein = p
                #max_score = 0
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #if hsp.score > max_score: max_score = hsp.score
                tmp_protein.list_of_shared_conserv_proteome_peps += Protein.Protein.hsp_match_parser(hsp.match,
                                                                                                     hsp.query,
                                                                                                     parsing_window_size=minlength, 
                                                                                                     max_sub=substitution, 
                                                                                                     max_mismatch=mismatch)
        if len(tmp_protein.list_of_shared_human_peps) == 0:
            outfile.write("\nNo shared peptides for " + query_name)   
        else:
            outfile.write("\nList of shared peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
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
