# #!/usr/bin/python3

import tmhmm                                                # to predict transmembrane domains
from Bio import SeqIO, pairwise2                            # to handle .fasta files and to align sequences
from Bio.Blast.Applications import NcbiblastpCommandline    # autoimmunity module, to query the local sapiens database
from Bio.ExPASy.ScanProsite import scan, read               # scan prosite functional domains information
from Bio.Blast import NCBIXML                               # to parse the peptides
from NERVE import Protein                                   # to contain proteins information
from NERVE.function import deep_fri, retrieve_entry_function
from NERVE.cello_scraper import *
from NERVE.quality_control import quality_control
from NERVE.proteome_downloader import proteome_downloader
import numpy as np                                          # to handle vectors and matrices (machine / deep learning modules)
import tensorflow                                           # for deep learning modules
from tensorflow import keras                                # to use the spaan model (to predict the probability of a protein to be an adhesin)
import pandas                                               # to read mhcpep.csv and write the final report
from spaan.data_processing import *                         # to extract proteins features for the spaan model
import sys
import urllib
import os
from shutil import rmtree
import logging

import argparse
from typing import NamedTuple

def dir_path(path):
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path

class Args(NamedTuple):
    '''Command-line arguments'''
    annotation:bool
    e_value:float
    gram:str
    minlength:int
    mismatch:int
    mouse:bool
    mouse_peptides_sum_limit:float
    proteome1:str
    proteome2:str
    p_ad_extracellular_filter:float
    p_ad_no_citoplasm_filter:float
    padlimit:float
    razor:bool
    razlen:int
    sapiens_peptides_sum_limit:float
    select:bool
    substitution:float
    transmemb_doms_limit:int
    verbose:int
    virlimit:float
    virulent:bool
    
    working_dir:str
    NERVE_dir:str
    
def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description="Run vaccine candidate prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-a',
                        metavar='--annotation', 
                        help="Activation or deactivation of annotation module, to retrieve info about protein functions. By default, this module is active. Type True or False to activate or deactivate it, respectively.",
                        type=bool,
                        required=False,
                        const=True
                        )
    parser.add_argument('-ev',
                        metavar='--e_value', 
                        help="Set expect-value used in blastp for immunity modules. For example: -ev=0.0001 or -e_value=0.0001",
                        type=float,
                        const=1e-10,
                        required=False,
                        )
    parser.add_argument('-g',
                        metavar='--gram', 
                        help="gram stain of the procaryotic pathogen of interest. Type p for positive or n for negative. For example: -g=p or -gram=p",
                        type=str,
                        required=True,
                        )
    parser.add_argument('-ml',
                        metavar='--minlength', 
                        help="minimal length required for shared peptides to be extracted in comparison analyses versus human and/or mouse. For example: -ml=10, -minlength=10",
                        type=int,
                        const=9,
                        required=False,
                        )
    parser.add_argument('-mm',
                        metavar='--mismatch', 
                        help="maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -mm=2, -mismatch=2",
                        type=int,
                        const=1,
                        required=False,
                        )
    parser.add_argument('-m',
                        metavar='--mouse', 
                        help="Activation or deactivation of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the eventual shared peptides is carried out as in the autoimmunity module. Type True or False to activate or deactivate it, respectively. For example: -m=True or -mouse=True",
                        type=bool,
                        const=False,
                        required=False,
                        )
    parser.add_argument('-mouse_peptides_sum_limit',
                        metavar='--mouse_peptides_sum_limit', 
                        help="Parameter calculated in mouse module and need to be used in select module. Protein with (sum of shared peptides of the i-protein with mouse proteins/number of aminoacids of the i-protein)<=0.15 and with absence of match mhc-I and Mhc-II mouse ligands are selected. It needs to be calulated for each protein",
                        type=float,
                        const=0.15,
                        required=False,
                        )
    parser.add_argument('-p1',
                        metavar='--proteome1', 
                        help='Path to proteome, alternative Uniprot proteome ID (https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-p2',
                        metavar='--proteome2', 
                        help='Path to proteome, alternative Uniprot proteome ID (https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=False,
                        )
    parser.add_argument('-p_ad_extracellular_filter',
                        metavar='--p_ad_extracellular_filter', 
                        help="Parameter of select module. Extracellular proteins with a pad lower than 0.38 are discarded",
                        type=float,
                        const=0.38,
                        required=False,
                        )
    parser.add_argument('-p_ad_no_citoplasm_filter',
                        metavar='--p_ad_no_citoplasm_filter', 
                        help="Parameter of select module. Non-cytoplasmic Proteins with a pad lower than 0.46 are discarded.",
                        type=float,
                        const=0.46,
                        required=False,
                        )    
    parser.add_argument('-padlimit',
                        metavar='--padlimit', 
                        help="Set the PAD value cut-off for proteins with 'Unknown' localization in the select module. Thus, these proteins with a PAD value < cut-off are discarded. Set a number between 0 and 1. For example: -pl=0.90, -padlimit=0.90",
                        type=float,
                        const=0.85,
                        required=False,
                        )
    parser.add_argument('-rz',
                        metavar='--razor', 
                        help="Activation or deactivation of the loop-razor module. This module allows you to recover protein vaccine candidates, with more than 2 transmembrane domains, that would be discarded in the last module. The longest loop with minimum 50 aa will replace the original protein sequence for next NERVE steps, if it is present. Type True or False to activate or deactivate it, respectively. For example: -rz=True or -razor=True",
                        type=bool,
                        const=False,
                        required=False,
                        )
    parser.add_argument('-rl',
                        metavar='--razlen', 
                        help="Set minimal length of loop considered in loop-razor module. For example: -rl=70 or -razlen=70",
                        type=int,
                        const=50,
                        required=False,
                        )
    parser.add_argument('-s',
                        metavar='--select', 
                        help="Activation or deactivation of select module, which filters PVC from proteome1. Type 'True' or 'False' to activate or deactivate it, respectively. For example: -s='False' or -select='False'",
                        type=bool,
                        const=True,
                        required=False,
                        )
    parser.add_argument('-ss',
                        metavar='--substitution', 
                        help="Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -ss=4, -substitution=4",
                        type=int,
                        const=3,
                        required=False,
                        )
    parser.add_argument('-transmemb_doms_limit',
                        metavar='--transmemb_doms_limit', 
                        help="Parameter of select module. Proiteins with trasmembrane domains>=3 are discarded.",
                        type=int,
                        const=3,
                        required=True,
                        )
    parser.add_argument('-v',
                        metavar='--verbose', 
                        help='',
                        type=int,
                        const=10,
                        required=False,
                        )
    parser.add_argument('-vl',
                        metavar='--virlimit', 
                        help="Fix a cut-off value for NERVirulent in the select module. Set a number between 0 and 1. (default = 0.5) For example: -vl=0.60 -virlimit=0.60",
                        type=float,
                        const=0.5,
                        required=True,
                        )    
    parser.add_argument('-virulent',
                        metavar='--virulent', 
                        help="Activation or deactivation of NERVirulent module, involved in the prediction of the probability of being a virulence factor through protein sequence analysis. Type True or False to activate or deactivate it, respectively. For example: -virulent=True",
                        type=bool,
                        const=False,
                        required=False,
                        )
     
    parser.add_argument('-working_dir',
                        metavar='--working_dir', 
                        help='Working directory',
                        type=dir_path,
                        required=False,
                        const='./'
                        )
    parser.add_argument('-NERVE_dir',
                        metavar='--NERVE_dir', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        const='../NERVE'
                        )
    
    args = parser.parse_args()
    return Args(args.annotation, args.e_value, args.gram, args.minlength, args.mismatch,
                args.mouse, args.mouse_peptides_sum_limit, args.proteome1, args.proteome2, 
                args.p_ad_extracellular_filter, args.p_ad_no_citoplasm_filter, args.padlimit, args.razor, 
                args.razlen, args.select, args.substitutionargs.transmemb_doms_limit, args.verbose, args.virlimit, 
                args.virulent, args.working_dir, args.NERVE_dir)

def main():
    args=get_args()
    # define log file
    logging.basicConfig(filename='.log',
                        filemode='w',
                        level=logging.DEBUG if args.verbose==10 else logging.CRITICAL,
                        force=True)
    # check input and download proteome:
    if os.path.isdir(os.path.join(working_dir, args.proteome1)) == False:
        try:
            proteome_downloader(args.proteome1, filename=os.path.join(args.working_dir,'proteome1.fasta'))
        except Exception as e:
            raise ValueError(f'{args.proteome1} rised the following error:\n{e}')
        logging.debug(f'{args.proteome1} succesfully downloaded')
    if proteome2:
        if os.path.isdir(os.path.join(working_dir, args.proteome2)) == False:
            try:
                proteome_downloader(args.proteome2, filename=os.path.join(args.working_dir,'proteome2.fasta'))
            except:
                raise logging.error(f'{args.proteome2} rised the following error:\n{e}')
            logging.debug(f'{args.proteome2} succesfully downloaded')
            
if __name__ == "__main__":
    main()
