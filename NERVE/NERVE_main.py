#!/usr/bin/python3
"""Run NERVE, reverse vaccinology software"""

#import sys
#sys.path.insert(0,'/NERVE')
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # disable most of the warnings 
import time
import requests
import argparse
from typing import NamedTuple
import tmhmm                                                # to predict transmembrane domains
from Bio import SeqIO, pairwise2                            # to handle .fasta files and to align sequences
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline    # autoimmunity module, to query the local sapiens database
from Bio.ExPASy.ScanProsite import scan, read               # scan prosite functional domains information
from Bio.Blast import NCBIXML                               # to parse the peptides
from NERVE import Protein                                   # to contain proteins information
from NERVE.function import deep_fri, retrieve_entry_function
from NERVE.cello_scraper import *
#from NERVE.proteome_downloader import proteome_downloader
import numpy as np                                          # to handle vectors and matrices (machine / deep learning modules)
import tensorflow                                           # for deep learning modules
from tensorflow import keras                                # to use the spaan model (to predict the probability of a protein to be an adhesin)
import pandas                                               # to read mhcpep.csv and write the final report
from spaan.data_processing import *                         # to extract proteins features for the spaan model
import sys
import urllib
from shutil import rmtree
import logging
import subprocess
import json
from operator import attrgetter

def dir_path(path:str)->str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path

class Args(NamedTuple):
    '''Command-line arguments'''
    annotation:str
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
    select:bool
    substitution:float
    subcell:str
    transmemb_doms_limit:int
    virlimit:float
    virulent:bool
    working_dir:str
    NERVE_dir:str
    iFeature_dir:str
    DeepFri_dir:str
    def print_args(self):
        return (f'''annotation: {self.annotation}, e_value: {self.e_value}, gram: {self.gram},
                minlength: {self.minlength}, mismatch: {self.mismatch}, mouse: {self.mouse},
                mouse_peptides_sum_limit: {self.mouse_peptides_sum_limit}, proteome1: {self.proteome1},
                proteome2: {self.proteome2}, p_ad_extracellular_filter: {self.p_ad_extracellular_filter},
                padlimit: {self.padlimit}, razor: {self.razor}, razlen: {self.razlen}, select: {self.select},
                substitution: {self.substitution}, subcell: {self.subcell}, transmemb_doms_limit: {self.transmemb_doms_limit},
                virlimit: {self.virlimit}, virulent: {self.virulent}, working_dir: {self.working_dir},
                NERVE_dir: {self.NERVE_dir}, iFeature_dir: {self.iFeature_dir},  DeepFri_dir: {self.DeepFri_dir}''')
    
def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description="Run vaccine candidate prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-a','--annotation',
                        metavar='\b', 
                        help="Activation or deactivation of annotation module, to retrieve info about protein functions. By default, this module is active. Type True or False to activate or deactivate it, respectively.",
                        type=str,
                        required=False,
                        default="True"
                        )
    parser.add_argument('-ev','--e_value',
                        metavar='\b', 
                        help="Set expect-value used in blastp for immunity modules. For example: -ev=0.0001 or -e_value=0.0001",
                        type=float,
                        default=1e-10,
                        required=False,
                        )
    parser.add_argument('-g','--gram',
                        metavar='\b', 
                        help="gram stain of the procaryotic pathogen of interest. Type p for positive or n for negative. For example: -g=p or -gram=p",
                        type=str,
                        required=True,
                        )
    parser.add_argument('-ml','--minlength',
                        metavar='\b', 
                        help="minimal length required for shared peptides to be extracted in comparison analyses versus human and/or mouse. For example: -ml=10, -minlength=10",
                        type=int,
                        default=9,
                        required=False,
                        )
    parser.add_argument('-mm','--mismatch',
                        metavar='\b', 
                        help="maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -mm=2, -mismatch=2",
                        type=int,
                        default=1,
                        required=False,
                        )
    parser.add_argument('-m','--mouse',
                        metavar='\b', 
                        help="Activation or deactivation of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the eventual shared peptides is carried out as in the autoimmunity module. Type True or False to activate or deactivate it, respectively. For example: -m=True or -mouse=True",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-mpsl','--mouse_peptides_sum_limit',
                        metavar='\b', 
                        help="Parameter calculated in mouse module and need to be used in select module. Protein with (sum of shared peptides of the i-protein with mouse proteins/number of aminoacids of the i-protein)<=0.15 and with absence of match mhc-I and Mhc-II mouse ligands are selected. It needs to be calulated for each protein",
                        type=float,
                        default=0.15,
                        required=False,
                        )
    parser.add_argument('-p1','--proteome1',
                        metavar='\b', 
                        help='Path to proteome, alternative Uniprot proteome ID (https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-p2','--proteome2',
                        metavar='\b', 
                        help='Path to proteome, alternative Uniprot proteome ID (https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=False,
                        )
    parser.add_argument('-paefilter','--p_ad_extracellular_filter',
                        metavar='\b', 
                        help="Parameter of select module. Extracellular proteins with a pad lower than 0.38 are discarded",
                        type=float,
                        default=0.38,
                        required=False,
                        )
    parser.add_argument('-pacfilter','--p_ad_no_citoplasm_filter',
                        metavar='\b', 
                        help="Parameter of select module. Non-cytoplasmic Proteins with a pad lower than 0.46 are discarded.",
                        type=float,
                        default=0.46,
                        required=False,
                        )    
    parser.add_argument('-pl','--padlimit',
                        metavar='\b',
                        help="Set the PAD value cut-off for proteins with 'Unknown' localization in the select module. Thus, these proteins with a PAD value < cut-off are discarded. Set a number between 0 and 1. For example: -pl=0.90, -padlimit=0.90",
                        type=float,
                        default=0.85,
                        required=False,
                        )
    parser.add_argument('-rz','--razor',
                        metavar='\b', 
                        help="Activation or deactivation of the loop-razor module. This module allows the recovery of protein vaccine candidates, with more than 2 transmembrane domains, that would otherwise be discarded in the last module. The longest loop with minimum 50 aa will replace the original protein sequence for following NERVE steps, if it is present. Type True or False to activate or deactivate it, respectively. For example: -rz=True or -razor=True",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-rl','--razlen',
                        metavar='\b', 
                        help="Set minimal length of loop considered in loop-razor module. For example: -rl=70 or -razlen=70",
                        type=int,
                        default=50,
                        required=False,
                        )
    parser.add_argument('-s','--select',
                        metavar='\b', 
                        help="Activation or deactivation of select module, which filters PVC from proteome1. Type 'True' or 'False' to activate or deactivate it, respectively. For example: -s='False' or -select='False'",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-ss','--substitution',
                        metavar='\b', 
                        help="Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -ss=4, -substitution=4",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-su','--subcell',
                        metavar='\b', 
                        help="Subcellular localization prediction program: 'cello' or 'psortb'",
                        type=str,
                        default="psortb",
                        required=False,
                        )
    parser.add_argument('-tdl','--transmemb_doms_limit',
                        metavar='\b', 
                        help="Parameter of select module. Proiteins with trasmembrane domains>=3 are discarded.",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-vl','--virlimit',
                        metavar='\b', 
                        help="Fix a cut-off value for NERVirulent in the select module. Set a number between 0 and 1. (default = 0.5) For example: -vl=0.60 -virlimit=0.60",
                        type=float,
                        default=0.5,
                        required=False,
                        )    
    parser.add_argument('-vir','--virulent',
                        metavar='\b', 
                        help="Activation or deactivation of NERVirulent module, involved in the prediction of the probability of being a virulence factor through protein sequence analysis. Type True or False to activate or deactivate it, respectively. For example: -virulent=True",
                        type=str,
                        default="True",
                        required=False,
                        )
     
    parser.add_argument('-wd','--working_dir',
                        metavar='\b', 
                        help='path to working directory. If not existing, a working directory with the given path is created',
                        type=str,
                        required=False,
                        default='./'
                        )
    parser.add_argument('-nd','--NERVE_dir',
                        metavar='\b', 
                        help='path to NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../NERVE'
                        )
    parser.add_argument('-id','--iFeature_dir',
                        metavar='\b', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../iFeature'
                        )
    parser.add_argument('-dfd','--DeepFri_dir',
                        metavar='\b', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../DeepFri'
                        )
    
    
    args = parser.parse_args()
    return Args(args.annotation, args.e_value, args.gram, args.minlength, args.mismatch,
                args.mouse, args.mouse_peptides_sum_limit, args.proteome1, args.proteome2, 
                args.p_ad_extracellular_filter, args.p_ad_no_citoplasm_filter, args.padlimit, args.razor, 
                args.razlen, args.select, args.substitution, args.subcell, args.transmemb_doms_limit, args.virlimit, 
                args.virulent, args.working_dir, args.NERVE_dir, args.iFeature_dir, args.DeepFri_dir)


def main():
    """Runs NERVE"""
    # record time:
    nerve_start = time.time()
    args=get_args()
    print("Start NERVE 1.5")
    
    # init workdir:
    if args.working_dir[-1] != '/':
        args=args._replace(working_dir=args.working_dir+'/')
    # create working directory if does not exist
    if os.path.isdir(args.working_dir)==False:
        os.makedirs(args.working_dir)
    # define log file
    logging.basicConfig(filename=os.path.join(args.working_dir, 'logfile.log'),
                        filemode='w',
                        level=logging.DEBUG,
                        force=True)
    logging.debug(f'Running NERVE with the following parameters:\n{args.print_args()}')    
    # check input and download proteome:
    if os.path.isfile(args.proteome1)==True:
        logging.debug(f'{args.proteome1} found as {args.proteome1}')
        # repath proteome as absolute path
        args = args._replace(proteome1=os.path.join('/', os.path.relpath(args.proteome1, start = '/')))
    elif os.path.isfile(os.path.join(args.working_dir, args.proteome1)) == True:
        logging.debug(f'{args.proteome1} was found in {args.working_dir}')
        args = args._replace(proteome1=os.path.join(args.working_dir, args.proteome1))
    else:
        logging.debug(f'{args.proteome1} is not a file, download from Uniprot.')
        try:
            proteome_downloader(args.working_dir, args.proteome1, filename=os.path.join(args.working_dir,\
                                                                                        'proteome1.fasta'))
        except Exception as e:
            raise ValueError(f'{args.proteome1} rised the following error:\n{e}')
        logging.debug(f'{args.proteome1} successfully downloaded')
        args = args._replace(proteome1=os.path.join(args.working_dir,'proteome1.fasta'))
    
    if args.proteome2:
        if os.path.isfile(args.proteome2)==True:
            logging.debug(f'{args.proteome2} found as {args.proteome2}')
            # repath proteome as absolute path
            args = args._replace(proteome2=os.path.join('/', os.path.relpath(args.proteome2, start = '/')))
        elif os.path.isfile(os.path.join(args.working_dir, args.proteome2)) == True:
            logging.debug(f'{args.proteome2} was found in {args.working_dir}')
            args = args._replace(proteome2=os.path.join(args.working_dir, args.proteome2))
        else:
            logging.debug(f'{args.proteome2} is not a file, download from Uniprot.')
            try:
                proteome_downloader(args.working_dir, args.proteome2, filename=os.path.join(args.working_dir,\
                                                                                            'proteome2.fasta'))
            except:
                raise logging.error(f'{args.proteome2} rised the following error:\n{e}')
            logging.debug(f'{args.proteome2} successfully downloaded')
            args = args._replace(proteome2=os.path.join(args.working_dir,'proteome2.fasta'))   
    print("10% done")
    
    # run quality control
    start=time.time()
    logging.debug(f'Start quality control of proteome1 ({args.proteome1})')
    # during the quality control, upload sequences from proteome1
    list_of_fasta_proteins, proteome1_new_path=quality_control(args.proteome1, args.working_dir, upload=True)
    # update input path of proteome1
    args=args._replace(proteome1=proteome1_new_path)
    logging.debug(f'Finish quality control of proteome1. Updated path: ({args.proteome1})')
    if args.proteome2:
        logging.debug(f'Start quality control of proteome2 ({args.proteome2})')
        proteome2_new_path=quality_control(args.proteome2, args.working_dir)
        # update proteome2 new path 
        args=args._replace(proteome2=proteome2_new_path)
        logging.debug(f'Finish quality control of proteome2. Updated path: ({args.proteome2})')
    logging.debug(f'Extract protein sequences and IDs from proteome1')
    list_of_proteins = []
    for p in list_of_fasta_proteins:
        p_id = str(p.name)
        p_seq = str(p.seq)
        list_of_proteins.append(Protein.Protein(p_id, p_seq))
    end=time.time()
    logging.debug(f'{len(list_of_fasta_proteins)} proteins loaded in {end-start} seconds')
            
    # subcellular localization prediction
    if args.subcell=='cello':
        start=time.time()
        logging.debug("Subcelloc start with cello...")
        list_of_proteins=cello(list_of_proteins, args.working_dir, args.gram, args.proteome1)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    
    if args.subcell=='psortb':
        start=time.time()
        logging.debug("Subcelloc start with psortb...")
        list_of_proteins=psortb(list_of_proteins, args.working_dir, args.gram, args.proteome1)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("20% done")
    
    # Adhesin
    logging.debug("Adhesin start...")
    start=time.time()
    list_of_proteins=adhesin(list_of_proteins, args.working_dir, args.NERVE_dir)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("30% done")
    
    # Tmhelices
    logging.debug("Tmhelices start...")
    start=time.time()
    list_of_proteins=Tmhelices(list_of_proteins, args.working_dir)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("40% done")
    
    # Razor
    if args.razor=="True":
        logging.debug("Loop-razor start...")
        start=time.time()
        list_of_proteins=razor(list_of_proteins, args.working_dir, args.transmemb_doms_limit, args.razlen)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("50% done")
    
    # Autoimmunity
    logging.debug("Autoimmunity start...")
    start=time.time()
    list_of_proteins=autoimmunity(list_of_proteins, args.proteome1, args.working_dir, args.NERVE_dir, args.e_value, 
                                  args.minlength, args.mismatch, args.substitution)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("60% done")
    
    # Mouse immunity
    if args.mouse=="True":
        start=time.time()
        logging.debug("Mouse immunity start...")
        list_of_proteins=mouse(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, args.proteome1,
                               args.minlength, args.substitution, args.mismatch)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    
    # Conservation
    #logging.debug(f'list of proteins before conservation:\n{list_of_proteins}')
    if args.proteome2:
        start=time.time()
        logging.debug("Conservation start...")
        list_of_proteins=conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, 
                                      args.proteome1, args.proteome2,
                                      args.minlength, args.substitution, args.mismatch)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("70% done")
        
    # Virulence
    #logging.debug(f'list of proteins before virulence:\n{list_of_proteins}')
    if args.virulent=="True":
        start=time.time()
        logging.debug("Virulence start...")
        list_of_proteins=virulence(list_of_proteins, args.working_dir, args.iFeature_dir, args.proteome1, args.NERVE_dir)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("80% done")
        
    # annotation
    #logging.debug(f'list of proteins before annotation:\n{list_of_proteins}')
    if args.annotation=="True":
        start=time.time()
        logging.debug("Annotation start...")
        list_of_proteins=annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("90% done")
    
    # select
    final_proteins=list_of_proteins
    if args.select=="True":
        logging.debug("Select start...")
        start=time.time()
        final_proteins=select(list_of_proteins, args.p_ad_no_citoplasm_filter, args.p_ad_extracellular_filter, 
               args.transmemb_doms_limit, args.padlimit, args.mouse, 
               args.mouse_peptides_sum_limit, args.virlimit, args.virulent)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))

    if args.virulent=="True":
        final_proteins.sort(key=lambda p: p.p_vir, reverse=True)
    
    # return .csv outputs
    output(final_proteins, os.path.join(args.working_dir, 'vaccine_candidates.csv'))
    # collect discarded proteins
    final_proteins_names = [p.id for p in final_proteins]
    discarded_proteins = [p for p in list_of_proteins if p.id not in final_proteins_names]
    output(discarded_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'))
    
    nerve_end = time.time()
    logging.debug("Done: NERVE has finished its analysis in: {:.4f} seconds".format(nerve_end-nerve_start))
    print("100% done")
    print("End NERVE computation successfully.")
    
def bashCmdMethod(bashCmd):
    """Run bash commands
    param: bashCmd: bash command to be run"""
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

class protein_element:
    """Class to handle fasta file elements similarly to the biopython fasta file parser"""
    def __init__(self, name, seq):
        self.name=name
        self.seq=seq

def proteome_downloader(working_dir, proteome_id, filename='input_proteome.fasta', output_dir=os.getcwd(), format_ = "fasta") -> None:
    """Downloads proteome from uniprot database into output multifasta file
    param: proteome_id: uniprot unique proteome id, not case-sensitive
    param: output_dir: output directory (default: current directory)
    param: format: uniprot API required format (default:fasta)
    param: filename: output proteome filename (default: input_proteome.fasta)
    """
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    try: 
        url=f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format={format_}&query=%28proteome%3A{proteome_id}%29'
        response = requests.get(url, stream = True)
        text_file = open(os.path.join(output_dir, filename), 'wb')
        for chunk in response.iter_content(chunk_size=1024):
              text_file.write(chunk)
        # raise an AssertionError if the given proteome ID is not valid
        assert text_file.tell() > 0, 'First download attempt failed'
        text_file.close()
    except AssertionError: 
        logging.debug(f'28proteome is not a valid building block')
        try:
            url=f'https://rest.uniprot.org/uniprotkb/stream?format={format_}&query=%28proteome%3A{proteome_id}%29'
            response = requests.get(url, stream = True)
            text_file = open(os.path.join(output_dir, filename), 'wb')
            for chunk in response.iter_content(chunk_size=1024):
                  text_file.write(chunk)
            # raise an AssertionError if the given proteome ID is not valid
            assert text_file.tell() > 0, 'Second download attempt failed'
            text_file.close()
        except AssertionError:
            logging.debug(f'avoiding compression is not a valid building block')
            try:
                url=f'https://rest.uniprot.org/uniparc/stream?format={format_}&query=%28upid%3A{proteome_id}%29'
                response = requests.get(url, stream = True)
                text_file = open(os.path.join(output_dir, filename), 'wb')
                for chunk in response.iter_content(chunk_size=1024):
                      text_file.write(chunk)
                # raise an AssertionError if the given proteome ID is not valid
                assert text_file.tell() > 0, 'Third download attempt failed'
                text_file.close()
            except AssertionError as e:
                logging.debug(f'28upid is not a valid building block')
                print(f'Unable to download proteome {proteome_id} due to invalid proteome ID Uniprot API failure or wrong path. In case of uniprot API failure provide the proteome as file.')
                raise SystemExit(e)
    return None

def proteome_uploader(infile:str)->list:
    """Function to read and parse fasta files. Bio SeqIO is not suitable because it chops sequence names. It 
    will be used only to validate fasta file format.
    param: infile: path to fasta file"""
    class protein_element:
        """Class to handle fasta file elements similarly to the biopython fasta file parser"""
        def __init__(self, name, seq):
            self.name=name
            self.seq=seq
    proteome_elements=[]
    proteome_data={}
    infile=open(infile, 'r').readlines()
    for i in range(len(infile)):
        name=infile[i].strip()[1:]
        if infile[i].startswith('>'):
            proteome_data[name]=''
        if infile[i].startswith('>')==False:
            proteome_data[name]+=infile[i].strip()
    for element in proteome_data:
        proteome_elements.append(protein_element(element, proteome_data[element]))
    return proteome_elements

def is_fasta(filename:str):
    """Function that rise an error if the format is not .fasta.
    param: filename: path to fasta file"""
    with open(filename, "r") as handle:
        fasta = list(SeqIO.parse(handle, "fasta"))
        # biopython silently fails if the format is not fasta returning an empty generator
        # any() returns False if the list is empty
        if any(fasta) == True:
            fasta=proteome_uploader(filename)
            return fasta
        else:
            raise ValueError(f'{filename} is not in fasta format')

def proteome_uploader(infile:str)->list:
    """Function to read and parse fasta files. Bio SeqIO is not suitable because it chops sequence names. It 
    will be used only to validate fasta file format.
    param: infile: path to fasta file"""
    class protein_element:
        """Class to handle fasta file elements similarly to the biopython fasta file parser"""
        def __init__(self, name, seq):
            self.name=name
            self.seq=seq
    proteome_elements=[]
    proteome_data={}
    infile=open(infile, 'r').readlines()
    for i in range(len(infile)):
        # use .strip() to remove initial spaces
        if infile[i].strip().startswith('>'):
            # upload protein name
            name=infile[i].strip()[1:]
            proteome_data[name]=''
        # upload sequence
        if infile[i].strip().startswith('>')==False and infile[i].strip()!='':
            # in case of wrong name input:
            #if name=="":
            #    raise ValueError('Encountered a sequence with wrong name formatting. Possible presence of spaces before ">" symbol.')
            proteome_data[name]+=infile[i].strip()
    for element in proteome_data:
        proteome_elements.append(protein_element(element, proteome_data[element]))
    return proteome_elements

def is_fasta(filename:str):
    """Function that rise an error if the format is not .fasta.
    param: filename: path to fasta file"""
    with open(filename, "r") as handle:
        fasta = list(SeqIO.parse(handle, "fasta"))
        # biopython silently fails if the format is not fasta returning an empty generator
        # any() returns False if the list is empty
        if any(fasta) == True:
            fasta=proteome_uploader(filename)
            return fasta
        else:
            raise ValueError(f'{filename} is not in fasta format')
            
def dir_path(path:str)->str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path

def quality_control(path_to_fasta:str, working_dir:str, upload=False)->dir_path:
    """
    Remove sequences with non-canonical aminoacid symbols. U (Se-Cys) is substituted with C (Cys). Returns
    {working_dir}/discarded_sequences_{input_fasta_file_name}.fasta containing discarded sequences and\
    {working_dir}/cleaned_{input_fasta_file_name}.fasta with cleaned sequences
    param: path_to_fasta: full path to fasta file containing the proteome with .fasta extension. Input fasta file\
    will not be overwritten
    param: working_dir: working directory, were cleaned fasta file will be saved
    param: upload: if True the function returns a list containing the filtered sequences as protein_element objects
    output: path to the cleaned fasta file named as {working_dir}/cleaned_{input_fasta_file_name}.fasta
    """
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    
    aa_dic = {'C': 'C', 'D': 'D', 'S': 'S', 'Q': 'Q', 'K': 'K', 'I': 'I', 'P': 'P', 'T': 'T', 'F': 'F', 'N': 'N', 
              'G': 'G', 'H': 'H', 'L': 'L', 'R': 'R', 'W': 'W', 'A': 'A', 'V': 'V', 'E': 'E', 'Y': 'Y', 'M': 'M', 
              'U':'C'}
    filtered_sequences, discarded_sequences = [],[]
    # control formatting
    fasta_list = is_fasta(path_to_fasta)    
    # filename needed to create the output file
    file_name = os.path.basename(path_to_fasta)
    output_file = os.path.join(working_dir, "_".join(["cleaned",file_name]))
    output_discarded_sequences=os.path.join(working_dir, "_".join(["discarded_sequences",file_name]))
    for record in fasta_list:
        flag = True
        new_seq = ''
        # check sequence
        for aa in str(record.seq):
            if aa not in aa_dic:
                flag = False
                logging.debug(f'Found non-canonical aminoacid "{aa}" in sequence {record.name}')
            elif aa == "U":
                logging.debug(f'Found non-canonical aminoacid "{aa}" (Selenocysteine) in sequence {record.name}, substituting to Cysteine')
                new_seq+=aa_dic[aa]
            else:
                new_seq+=aa_dic[aa]
        # check name
        if ">" in record.name:
            logging.debug(f'Found non-canonical character ">" in sequence name:\n{record.name}\nSubstituting with "*"')
            record.name=record.name.replace(">","*")
        record.seq=Seq(new_seq)
        
        if flag==True:
            filtered_sequences.append(record)
        else:    
            discarded_sequences.append(record)
            logging.debug(f'Sequence {record.name} has been discarded for the presence of non-canonical aminoacids.')     
    # output filtered overwriting input fasta file
    filename = open(output_file, 'w')
    for sequence in filtered_sequences:
        filename.write(f'>{str(sequence.name)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    #SeqIO.write(filtered_sequences, filename, "fasta")
    filename.close()
    # output discarded sequences
    filename = open(output_discarded_sequences, 'w')
    for sequence in discarded_sequences:
        filename.write(f'>{str(sequence.name)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    #SeqIO.write(discarded_sequences, filename, "fasta")
    filename.close()
    if upload==True:
        # return filtered sequnces and the new path
        return filtered_sequences, output_file
    # return the path to the cleaned file
    return output_file

def cello(list_of_proteins, working_dir, gram, proteome1)->list:
    "Run cello scraper for subcellular localization prediction"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    cello_text_output = cello_scraper(gram, proteome1)
    df = cello_output_parser(cello_text_output)
    df.columns = ["name", "prediction"]
    treshold = 1
    for index, row in df.iterrows():
        for p in list_of_proteins:
            if p.id == row["name"] and row["prediction"][0].reliability >= treshold:
                p.localization=row["prediction"]
                #p.localization = [pred.localization for pred in row["prediction"]]+[pred.reliability for pred in row["prediction"]]
    # save cello raw predictions
    #df.to_csv(os.path.join(working_dir, 'cello_predictions.csv'))
    return list_of_proteins

def psortb(list_of_proteins, working_dir, gram, proteome1)->list:
    """Psortb request sender and parser"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    class Localization:
        """class to store and handle localizations"""
        def __init__(self, localization, reliability):
            self.localization=str(localization)
            self.reliability=float(reliability)
    # read proteome1 file and pass it to json request
    infile=open(proteome1, 'r')
    proteome1="".join(infile)
    body={'gram':gram,'seq':proteome1}
    with open(os.path.join(working_dir, 'payload.json'), 'w') as f:
        json.dump(body, f)
    url="http://psortb:8080/execute"
    req = urllib.request.Request(url)
    req.add_header('Content-Type', 'application/json; charset=utf-8')
    jsondata = json.dumps(body)
    jsondataasbytes = jsondata.encode('utf-8')   # needs to be bytes
    req.add_header('Content-Length', len(jsondataasbytes))
    logging.debug("Sending request to psortb")
    response = urllib.request.urlopen(req, jsondataasbytes)
    data_json = json.loads(response.read())
    logging.debug(f'Psortb stdout:\n{data_json["stdout"]}\nPsortb stderr:\n{data_json["stderr"]}')
    logging.debug('Parsing psortb output')
    for entry in data_json['result'].split('-------------------------------------------------------------------------------\n\n'):
        split=entry.split('\n')
        id_=split[0][split[0].find('SeqID: ')+len('SeqID: '):].strip() # protein id
        # extract sublocalization predictions:
        if 'Localization Scores:' in entry:
            # get part of the output that define the predictions
            predictions=entry[entry.find('Localization Scores:')+len('Localization Scores:'):entry.find('Final')]
            # get list containing each prediction and its score
            predictions=(predictions.strip().split('\n'))
            localizations=[Localization(element.split()[0], element.split()[1]) for element in predictions]
            # sort localizations
            localizations = sorted(localizations, key=attrgetter('reliability'), reverse=True)
            # set unknown prediction
            if localizations[0].reliability <= 3.:
                localizations = [Localization('Unknown', 0)]
            for p in list_of_proteins:
                if p.id == id_:
                    p.localization=localizations
    # save psortb raw predictions
    #df.to_csv(os.path.join(working_dir, 'psortb_predictions.csv'))
    return list_of_proteins

def adhesin(list_of_proteins, working_dir, NERVE_dir)->list:
    "Runs adhesin predictions"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    # take the model from nerve but the methods from spaan, cause the model in spaan could be modified and tested
    model = keras.models.load_model(os.path.join(os.path.join(NERVE_dir, "models"), 'espaan_model.h5')) 
    for p in list_of_proteins:
        p.p_ad = float(model.predict([
                     np.array([aminoacids_frequencies(p.sequence)]),
                     np.array([multiplet_frequencies(p.sequence, 3)]),
                     np.array([multiplet_frequencies(p.sequence, 4)]),
                     np.array([multiplet_frequencies(p.sequence, 5)]),
                     np.array([dipeptide_frequencies(p.sequence)]),
                     np.array([charge_composition(p.sequence)]),
                     np.array([hydrophobic_composition(p.sequence)])
                 ]))
    return list_of_proteins

def Tmhelices(list_of_proteins, working_dir)->list:
    "Runs Tmhelices"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    for p in list_of_proteins:
        annotation, _ = tmhmm.predict(p.sequence)
        p.tmhmm_seq = annotation
        transmembrane_domains = 0
        for i in range(len(annotation)-1):
            if (annotation[i] == 'i' or annotation[i] == 'o') and annotation[i+1] == 'M':
                transmembrane_domains += 1
        p.transmembrane_doms = transmembrane_domains
    return list_of_proteins
        

def razor(list_of_proteins, working_dir, transmemb_doms_limit, razlen)->list:
    "Runs razor module"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    #logging.debug("Warning: razor uses X as an exclusive symbol to split the final protein. Check if X is used inside the protein sequence!")
    for p in list_of_proteins:
        if p.transmembrane_doms >= transmemb_doms_limit:
            longest_loop = max(p.provide_raw_loops(), key = lambda k: len(k))
            if len(longest_loop) > razlen:
                logging.debug(f"Substituting {str(p.id)} sequence with its longest loop.")
                p.original_sequence_if_razor = p.sequence
                p.sequence = longest_loop
    return list_of_proteins

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
    outfile=open(os.path.join(working_dir, 'autoimmunity_raw_output.txt'), 'w')
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
        if len(tmp_protein.list_of_shared_human_peps) == 0:
            outfile.write("\nNo immunogenic peptides for " + query_name)   
        else:
            outfile.write("\nList of immunogenic peptides for " + query_name + ": " +\
                          str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    outfile.close()
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
    mhcpep = pandas.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_sapiens.csv"), skipinitialspace=True)
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
    outfile=open(os.path.join(working_dir, 'mouse_immunity_raw_output.txt'), 'w')
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
        if len(tmp_protein.list_of_shared_human_peps) == 0:
            outfile.write("\nNo immunogenic peptides for " + query_name)   
        else:
            outfile.write("\nList of immunogenic peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    outfile.close()
    os.remove(os.path.join(working_dir, "mouse.xml")) # delete after the computation
    
    # store peptides from comparison with mouse recognized bacterial mhcpep
    mhcpep = pandas.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_mouse.csv"), skipinitialspace=True)
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
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(working_dir, 'compare_proteome/compare_proteome'), evalue=e_value, outfmt=5, out=os.path.join(working_dir,"comparison.xml")) # 5 is for xml 
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
        rmtree(os.path.join(working_dir, "compare_proteome"))
    return list_of_proteins

def virulence(list_of_proteins, working_dir, iFeature_dir, proteome1, NERVE_dir)->list:
    """Run virulent factor predictions"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type AAC --out {os.path.join(working_dir,'aac.out')}")
    bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type DPC --out {os.path.join(working_dir,'dpc.out')}")
    bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type CTDC --out {os.path.join(working_dir,'ctdc.out')}")
    bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type CTDT --out {os.path.join(working_dir,'ctdt.out')}")
    bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type CTDD --out {os.path.join(working_dir,'ctdd.out')}")
    extension = ".out"
    files = ["aac", "dpc", "ctdc", "ctdt", "ctdd"]
    datasets = [[] for el in files]
    for i in range(len(files)):
        with open(os.path.join(working_dir, files[i]+extension)) as f:
            lines = f.readlines()[1:]
            check_prot = 0
            for line in lines:
                information = line.split('\t')
                if not information[0] == list_of_proteins[check_prot].id:
                    logging.debug("Error in protein order! Return")
                datasets[i].append(np.array([float(el) for el in information[1:]]))
                check_prot += 1
        datasets[i] = np.array(datasets[i])
    labels = np.array([0. for _ in range(len(datasets[0]))])
    virulent_model = tensorflow.keras.models.load_model(os.path.join(os.path.join(NERVE_dir,"models"),\
                                                                     'virulent_classification_model.h5'))
    for i in range(len(datasets)):
        for j in range(len(datasets[i])):
            datasets[i][j] = np.array(datasets[i][j])
        datasets[i] = np.array(datasets[i])
    prediction = virulent_model.predict(datasets)
    for i in range(len(prediction)):
        list_of_proteins[i].p_vir = prediction[i][0]
    for file in files:
        os.remove(os.path.join(working_dir, file+extension)) # delete after the computation
    return list_of_proteins

def annotation(list_of_proteins, proteome1, working_dir, DeepFri_dir)->list:
    """Run protein function prediction"""
    deepfri_df = deep_fri(proteome1, DeepFri_dir, working_dir)
    for p in list_of_proteins:
        for index, row in deepfri_df.iterrows():
            if row['Protein'] in p.id:
                p.annotations = row['Function']
    return list_of_proteins

def select(list_of_proteins, p_ad_no_citoplasm_filter, p_ad_extracellular_filter, transmemb_doms_limit,
           padlimit, mouse, mouse_peptides_sum_limit, virlimit, virulent)->list:
    """Selects suitable candidate proteins for vaccine production"""
    final_list = []
    for protein in list_of_proteins:
        if protein.localization[0].localization == "Cytoplasmic": continue 
        #if protein.localization[0].reliability <= 3: continue
        if protein.p_ad < p_ad_no_citoplasm_filter and not protein.localization[0].localization == "Extracellular": continue 
        if protein.p_ad < p_ad_extracellular_filter and protein.localization[0].localization == "Extracellular": continue 
        if (protein.transmembrane_doms >= transmemb_doms_limit) and (protein.original_sequence_if_razor is None): continue
        if protein.sapiens_peptides_sum > .15: continue
        if len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens) >= 1: continue
        if (float(protein.localization[0].reliability) < 7.49) and (protein.p_ad < padlimit): continue
        #if (protein.localization[0].localization == "Unknown") and (protein.p_ad < padlimit): continue
        if mouse==True:
            if protein.mouse_peptides_sum > mouse_peptides_sum_limit: continue 
            if len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse) >= 1: continue 
        if virulent==True:
            if protein.p_vir < virlimit: continue
        final_list.append(protein)
    return final_list

def output(list_of_proteins, outfile):
    pd.DataFrame([[str(protein.id),
                 str("".join([str(protein.accession) if protein.accession!=None else ""])),
                 str(protein.length),
                 str(protein.transmembrane_doms),
                 str(protein.localization[0].localization),
                 str(protein.localization[0].reliability),
                 #str(", ".join([str(element) for element in protein.localization])),
                 str("".join([str(round(protein.p_vir,4)) if protein.p_vir!=None else ""])),
                 str("".join([str(round(protein.p_ad, 4)) if protein.p_ad!=None else ""])),
                 str("".join([str(round(protein.conservation_score, 4)) if protein.conservation_score!=None else ""])),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_human_peps if len(protein.list_of_shared_human_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_mouse_peps if len(protein.list_of_shared_mouse_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_conserv_proteome_peps if len(protein.list_of_shared_conserv_proteome_peps)>0])))),
                 str("".join([str(round(protein.sapiens_peptides_sum,4)) if protein.sapiens_peptides_sum!=None else "0"])),
                 str("".join([str(round(protein.mouse_peptides_sum,4)) if protein.mouse_peptides_sum!=None else "0"])),
                 str("".join([str(protein.annotations) if protein.annotations!=None else ""])),
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)))), 
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)))),  
                 str(protein.sequence),
                 str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor!=None else ""])),
                 str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""]))
                 ] for protein in list_of_proteins
                ], 
                columns= ['id ',
                    'uniprot_accession_code',
                    'length',
                    'transmembrane_doms',
                    'localization',
                    'localization score',
                    'virulence_probability',
                    'adhesin_probability',
                    'conservation_score',
                    'list_of_shared_human_peps',
                    'list_of_shared_mouse_peps',
                    'list_of_shared_conserv_proteome_peps',
                    'human_peptides_sum',
                    'mouse_peptides_sum',
                    'annotations',
                    'list_of_peptides_from_comparison_with_mhcpep_sapiens',
                    'list_of_peptides_from_comparison_with_mhcpep_mouse',
                    'sequence',
                    'original_sequence_if_razor',
                    'tmhmm_seq'
                     ]
                ).to_csv(outfile) 
    
if __name__ == "__main__":
    main()
