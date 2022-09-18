# #!/usr/bin/python3

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
    
    # set directories
    #current_dir = os.getcwd()
    #DeepFri_dir = os.path.join(current_dir, "DeepFri")
    
    # interrogate uniprot for function
    fastas = open(args.path_to_fastas, 'r')
    fasta_list = list(SeqIO.parse(fastas, "fasta"))
    outfile = os.path.join(working_dir, "deep_fri_input.fa")
    UniProt_df = retrieve_entry_function(fastas, outfile)
    
    DeepFri_df = deep_fri("deep_fri_input.fa", DeepFri_dir)
    
    # merge the two dfs
    final_df = pd.concat([UniProt_df, DeepFri_df])
    
    # add excluded proteins
    catched_proteins = final_df['Protein'].unique()
    listout = [[protein.id, 'Uknown function'] for protein in fasta_list if protein.id not in catched_proteins]
    missing_df = pd.DataFrame(listout, columns=['Protein', 'Function'])
    final_df = pd.concat([final_df, missing_df])
    final_df.to_csv(os.path.join(current_dir, 'annotation.csv'), index = False)
    #print(f'Results are saved as annotation.csv. This is a preview:\n{final_df}')

def deep_fri(path_to_fasta, DeepFri_dir, working_dir):
    current_dir = os.getcwd()
    #working_dir+="/" if working_dir[-1] != "/" else ""
    if working_dir[0] == ".":
        working_dir=working_dir[working_dir.find("/")+1:]
    # launch DeepFri
    os.chdir(DeepFri_dir)
    bashCmd = f'python3 predict.py --fasta_fn {path_to_fasta} -ont mf -o DeepFri'#{os.path.join(current_dir, working_dir)}'
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.chdir(current_dir)
    
    # parse DeepFri results
    DeepFri_df = DeepFriParser(os.path.join(working_dir, 'DeepFri_MF_predictions.csv'))
    return DeepFri_df

def retrieve_entry_function(fastas: str, outfile: str) -> pd.DataFrame: 
    '''Retrieves Uniprot annotated entry function of provided ACs and returns the data in a dictionary format
    param: fasta_list: list of Biopython Seq.IO objects containing protein ids and sequences to be evaluated;
    param: outfile: fasta file output containing sequences whose function was not found in the Uniprot
    output: output_df: data frame with Protein and Function fields where protein is the input fasta id and function is the annotated protein function from the uniprot database'''
    fasta_list = list(SeqIO.parse(fastas, "fasta"))
    dic, tdp, AC, = [{}, [], str()]
    # extract functions from uniprot
    for p in fasta_list:
        try:
            Id = p.id
            AC = p.id.split("|")[1]
            handle = urllib.request.urlopen(f"http://www.uniprot.org/uniprot/{AC}.xml")
            record = SeqIO.read(handle, "uniprot-xml")
            dic[Id] = [f'Uniprot annotations: {" | ".join([i for i in record.annotations["comment_function"]])}']
        except KeyError as e:
            if e.args[0] == 'comment_function':
                tdp.append(p)  
    # write outfile with proteins whose function was not found
    filename = open(outfile, 'w')
    SeqIO.write(tdp, filename, "fasta")
    filename.close()
    # create and rearrange dataframe
    output_df = pd.DataFrame.from_dict(dic, orient='index')
    if len(output_df.columns) == 2:
       output_df = output_df.reset_index()
       output_df.columns = ['Protein', 'Function']
       return output_df
    else:
       print("No Uniprot matches found")
       return pd.DataFrame([['', '']], columns=['Protein', 'Function'])

def DeepFriParser(path_to_infile: open) -> pd.DataFrame:
    '''Parses DeepFri input file
    param: path_to_infile: path to .csv output file produced by DeepFri program
    output: output_df: dataframe with Protein and Function fields where protein is the input protein fasta id and function is the | separated list of GO terms in descending order based on their score'''
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
        GOs = f'DeepFri predictions: {" | ".join([element.GO for element in sorted_row_tuples if element.score >=0.3])}'
        listout.append([Protein, GOs])
    output_df = pd.DataFrame(listout, columns = ['Protein', 'Function'])
    return output_df

if __name__ == "__main__":
    main()
