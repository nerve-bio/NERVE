#!/usr/bin/python3
"""Runs protein function prediction with DeepFri"""

import os, urllib
import pandas as pd
from Bio import SeqIO
from typing import NamedTuple
from operator import attrgetter
from code.Utils import run_command, build_id_index, match_protein

class Args(NamedTuple):
    '''Command-line arguments'''
    path_to_fastas:str
    DeepFri_dir:str
    working_dir:str

class DeepFriEntry(NamedTuple):
    '''Handles DeepFri entry'''
    id_:str
    score:float
    GO:str

def deep_fri(path_to_fasta, DeepFri_dir, working_dir):
    """Runs DeepFri prediction of protein function
       param: path_to_fasta: path to input fasta sequence file
       param: DeepFri_dir: path to /DeepFri/ directory
       param: working_dir: path to working directory
       output: pandas DataFrame containing fasta ID associated to protein function"""
    current_dir = os.getcwd()
    #working_dir+="/" if working_dir[-1] != "/" else ""
    if working_dir[0] == ".":
        working_dir=working_dir[working_dir.find("/")+1:]
    # launch DeepFri; a stale .csv would be parsed as this run's output
    output_path = os.path.join(working_dir, '_MF_predictions.csv')
    if os.path.exists(output_path):
        os.remove(output_path)
    run_command(['python3', 'predict.py', '--fasta_fn', path_to_fasta, '-ont', 'mf',
                 '-o', os.path.join(current_dir, working_dir)], cwd=DeepFri_dir)
    if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
        raise RuntimeError(f'DeepFri did not produce {output_path}')
    
    # parse DeepFri results
    DeepFri_df = DeepFriParser(output_path)
    return DeepFri_df

def DeepFriParser(path_to_infile: open, treshold=0.5) -> pd.DataFrame:
    '''Parses DeepFri input file
    param: path_to_infile: path to .csv output file produced by DeepFri program
    output: output_df: dataframe with Protein and Function fields where protein is the input protein fasta 
            id and function is the | separated list of GO terms in descending order based on their score. 
            Only elements with a score > treshold are considered valid
    treshold: float, DeepFri treshold (0.3)       
            '''
    # import the .csv file in a table-like format using pandas module 
    df = pd.read_csv(path_to_infile, comment='#')
    # rearrange information output
    listout = []
    for protein in df['Protein'].unique():
        subdf = df[df['Protein'] == protein]
        # get entries from every row in the form of named tuples
        row_tuples = [DeepFriEntry(id_=row["Protein"], score=row["Score"], GO=row["GO_term/EC_number name"]) for index, row in subdf.iterrows()]
        # sort named tuples based on score value in descending order
        Protein = row_tuples[0].id_
        sorted_row_tuples = sorted(row_tuples, key=attrgetter('score'), reverse=True)
        # get only gene ontology values
        GOs = f'{" | ".join([element.GO for element in sorted_row_tuples if element.score > treshold])}'
        listout.append([Protein, GOs])
    output_df = pd.DataFrame(listout, columns = ['Protein', 'Function'])
    return output_df
    
def annotation(list_of_proteins, proteome1, working_dir, DeepFri_dir)->list:
    """Run protein function prediction"""
    deepfri_df = deep_fri(proteome1, DeepFri_dir, working_dir)
    protein_index = build_id_index(list_of_proteins)
    for index, row in deepfri_df.iterrows():
        protein = match_protein(row['Protein'], protein_index, list_of_proteins)
        if protein is not None:
            protein.annotations = row['Function']
    for file in ['_MF_predictions.csv', '_MF_pred_scores.json']:
        os.remove(os.path.join(working_dir, file))
    return list_of_proteins
