#!/usr/bin/python3

"""Runs protein function prediction with DeepFri"""

import argparse, os, urllib, subprocess
import pandas as pd
from Bio import SeqIO
from typing import NamedTuple
from operator import attrgetter

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

def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description='Extract protein information from uniprot where accessible or predict with DeepFri',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-path_to_fastas',
                        metavar='--path-to-fastas', 
                        help='Path where sequences to evaluate are stored',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-DeepFri_dir',
                        metavar='--DeepFri_dir', 
                        help='Path where DeepFri repository is',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-working_dir',
                        metavar='--working_dir', 
                        help='Working directory path',
                        type=str,
                        required=True,
                        )
                        
    args = parser.parse_args()
    return Args(args.path_to_fastas, args.DeepFri_dir, args.working_dir)

def main() -> None:
    '''Executes functions of this module'''
    # get command-line arguments
    args = get_args()
    
    DeepFri_df = deep_fri(args.path_to_fasta, args.DeepFri_dir, args.working_dir)
    
    # add excluded proteins
    catched_proteins = DeepFri_df['Protein'].unique()
    listout = [[protein.id, 'Uknown function'] for protein in fasta_list if protein.id not in catched_proteins]
    missing_df = pd.DataFrame(listout, columns=['Protein', 'Function'])
    final_df = pd.concat([DeepFri_df, missing_df])
    final_df.to_csv(os.path.join(current_dir, 'annotation.csv'), index = False)
    print(f'Results are saved as annotation.csv. This is a preview:\n{final_df}')

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
    # launch DeepFri
    os.chdir(DeepFri_dir)
    bashCmd = f'python3 predict.py --fasta_fn {path_to_fasta} -ont mf -o {os.path.join(current_dir, working_dir)}'
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.chdir(current_dir)
    
    # parse DeepFri results
    DeepFri_df = DeepFriParser(os.path.join(working_dir, '_MF_predictions.csv'))
    return DeepFri_df

def DeepFriParser(path_to_infile: open) -> pd.DataFrame:
    '''Parses DeepFri input file
    param: path_to_infile: path to .csv output file produced by DeepFri program
    output: output_df: dataframe with Protein and Function fields where protein is the input protein fasta id and function is the | separated list of GO terms in descending order based on their score. Only elements with a score >0.5 are considered valid'''
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
        GOs = f'DeepFri predictions: {" | ".join([element.GO for element in sorted_row_tuples if element.score >0.5])}'
        listout.append([Protein, GOs])
    output_df = pd.DataFrame(listout, columns = ['Protein', 'Function'])
    return output_df

if __name__ == "__main__":
    main()
