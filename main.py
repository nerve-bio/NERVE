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
import subprocess

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
    select:bool
    substitution:float
    transmemb_doms_limit:int
    virlimit:float
    virulent:bool
    working_dir:str
    NERVE_dir:str
    iFeature_dir:str
    DeepFri_dir:str
    
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
                        default=True
                        )
    parser.add_argument('-ev',
                        metavar='--e_value', 
                        help="Set expect-value used in blastp for immunity modules. For example: -ev=0.0001 or -e_value=0.0001",
                        type=float,
                        default=1e-10,
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
                        default=9,
                        required=False,
                        )
    parser.add_argument('-mm',
                        metavar='--mismatch', 
                        help="maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -mm=2, -mismatch=2",
                        type=int,
                        default=1,
                        required=False,
                        )
    parser.add_argument('-m',
                        metavar='--mouse', 
                        help="Activation or deactivation of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the eventual shared peptides is carried out as in the autoimmunity module. Type True or False to activate or deactivate it, respectively. For example: -m=True or -mouse=True",
                        type=bool,
                        default=False,
                        required=False,
                        )
    parser.add_argument('-mpsl',
                        metavar='--mouse_peptides_sum_limit', 
                        help="Parameter calculated in mouse module and need to be used in select module. Protein with (sum of shared peptides of the i-protein with mouse proteins/number of aminoacids of the i-protein)<=0.15 and with absence of match mhc-I and Mhc-II mouse ligands are selected. It needs to be calulated for each protein",
                        type=float,
                        default=0.15,
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
    parser.add_argument('-paefilter',
                        metavar='--p_ad_extracellular_filter', 
                        help="Parameter of select module. Extracellular proteins with a pad lower than 0.38 are discarded",
                        type=float,
                        default=0.38,
                        required=False,
                        )
    parser.add_argument('-pacfilter',
                        metavar='--p_ad_no_citoplasm_filter', 
                        help="Parameter of select module. Non-cytoplasmic Proteins with a pad lower than 0.46 are discarded.",
                        type=float,
                        default=0.46,
                        required=False,
                        )    
    parser.add_argument('-pl',
                        metavar='--padlimit', 
                        help="Set the PAD value cut-off for proteins with 'Unknown' localization in the select module. Thus, these proteins with a PAD value < cut-off are discarded. Set a number between 0 and 1. For example: -pl=0.90, -padlimit=0.90",
                        type=float,
                        default=0.85,
                        required=False,
                        )
    parser.add_argument('-rz',
                        metavar='--razor', 
                        help="Activation or deactivation of the loop-razor module. This module allows you to recover protein vaccine candidates, with more than 2 transmembrane domains, that would be discarded in the last module. The longest loop with minimum 50 aa will replace the original protein sequence for next NERVE steps, if it is present. Type True or False to activate or deactivate it, respectively. For example: -rz=True or -razor=True",
                        type=bool,
                        default=False,
                        required=False,
                        )
    parser.add_argument('-rl',
                        metavar='--razlen', 
                        help="Set minimal length of loop considered in loop-razor module. For example: -rl=70 or -razlen=70",
                        type=int,
                        default=50,
                        required=False,
                        )
    parser.add_argument('-s',
                        metavar='--select', 
                        help="Activation or deactivation of select module, which filters PVC from proteome1. Type 'True' or 'False' to activate or deactivate it, respectively. For example: -s='False' or -select='False'",
                        type=bool,
                        default=True,
                        required=False,
                        )
    parser.add_argument('-ss',
                        metavar='--substitution', 
                        help="Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules. For example: -ss=4, -substitution=4",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-tdl',
                        metavar='--transmemb_doms_limit', 
                        help="Parameter of select module. Proiteins with trasmembrane domains>=3 are discarded.",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-vl',
                        metavar='--virlimit', 
                        help="Fix a cut-off value for NERVirulent in the select module. Set a number between 0 and 1. (default = 0.5) For example: -vl=0.60 -virlimit=0.60",
                        type=float,
                        default=0.5,
                        required=False,
                        )    
    parser.add_argument('-vir',
                        metavar='--virulent', 
                        help="Activation or deactivation of NERVirulent module, involved in the prediction of the probability of being a virulence factor through protein sequence analysis. Type True or False to activate or deactivate it, respectively. For example: -virulent=True",
                        type=bool,
                        default=False,
                        required=False,
                        )
     
    parser.add_argument('-wd',
                        metavar='--working_dir', 
                        help='Working directory',
                        type=dir_path,
                        required=False,
                        default='./'
                        )
    parser.add_argument('-nd',
                        metavar='--NERVE_dir', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../NERVE'
                        )
    parser.add_argument('-id',
                        metavar='--iFeature_dir', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../iFeature'
                        )
    parser.add_argument('-dfd',
                        metavar='--DeepFri_dir', 
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='../DeepFri'
                        )
    
    
    args = parser.parse_args()
    return Args(args.a, args.ev, args.g, args.ml, args.mm,
                args.m, args.mpsl, args.p1, args.p2, 
                args.paefilter, args.pacfilter, args.pl, args.rz, 
                args.rl, args.s, args.ss, args.tdl, args.vl, 
                args.vir, args.wd, args.nd, args.id, args.dfd)

def main():
    """Runs NERVE"""
    args=get_args()
    # define log file
    logging.basicConfig(filename=os.path.join(args.working_dir, 'logfile.log'),
                        filemode='w',
                        level=logging.DEBUG,
                        force=True)
    # check input and download proteome:
    logging.debug(f'Looking for {args.proteome1} in {args.working_dir}')
    if os.path.isfile(os.path.join(args.working_dir, args.proteome1)) == False:
        logging.debug(f'{args.proteome1} is not a file, download from Uniprot.')
        try:
            proteome_downloader(args.proteome1, filename=os.path.join(args.working_dir,'proteome1.fasta'))
        except Exception as e:
            raise ValueError(f'{args.proteome1} rised the following error:\n{e}')
        logging.debug(f'{args.proteome1} succesfully downloaded')
        args = args._replace(proteome1=os.path.join(args.working_dir,'proteome1.fasta'))
    else:
        logging.debug(f'{args.proteome1} was succesfully found in {args.working_dir}')
        args = args._replace(proteome1=os.path.join(args.working_dir, args.proteome1))
    if args.proteome2:
        logging.debug(f'Looking for {args.proteome2} in {args.working_dir}')
        if os.path.isfile(os.path.join(args.working_dir, args.proteome2)) == False:
            logging.debug(f'{args.proteome2} is not a file, download from Uniprot.')
            try:
                proteome_downloader(args.proteome2, filename=os.path.join(args.working_dir,'proteome2.fasta'))
            except:
                raise logging.error(f'{args.proteome2} rised the following error:\n{e}')
            logging.debug(f'{args.proteome2} succesfully downloaded')
            args = args._replace(proteome2=os.path.join(args.working_dir,'proteome2.fasta'))
        else:
            logging.debug(f'{args.proteome2} was succesfully found in {args.working_dir}')
            args = args._replace(proteome2=os.path.join(args.working_dir, args.proteome2))
    
    # run quality control
    logging.debug(f'Start quality control of proteome1 ({args.proteome1})')
    quality_control(args.proteome1, args.working_dir)
    logging.debug(f'Finish quality control of proteome1 ({args.proteome1})')
    if args.proteome2:
        logging.debug(f'Start quality control of proteome2 ({args.proteome2})')
        quality_control(args.proteome2, args.working_dir)
        logging.debug(f'Finish quality control of proteome2 ({args.proteome2})')
    
    # extract protein sequences and IDs:
    logging.debug(f'Extract protein sequences and IDs from proteome1')
    list_of_fasta_proteins = list(SeqIO.parse(args.proteome1, "fasta"))
    list_of_proteins = []
    for p in list_of_fasta_proteins:
        p_id = p.id
        p_seq = p.seq
        list_of_proteins.append(Protein.Protein(p_id, p_seq))
            
    # subcellular localization prediction
    logging.debug("Subcelloc start...")
    list_of_proteins=cello(list_of_proteins, args.working_dir, args.gram, args.proteome1)
    logging.debug("Done.")
    
    # Adhesin
    logging.debug("Adhesin start...")
    list_of_proteins=adhesin(list_of_proteins, args.working_dir, args.NERVE_dir)
    logging.debug("Done.")
    
    # Tmhelices
    logging.debug("Tmhelices start...")
    list_of_proteins=Tmhelices(list_of_proteins, args.working_dir)
    logging.debug("Done.")
    
    # Razor
    if args.razor:
        logging.debug("Loop-razor start...")
        list_of_proteins=razor(list_of_proteins, args.working_dir, args.transmemb_doms_limit, args.razor_len)
        logging.debug("Done.")
    
    # Autoimmunity
    logging.debug("Autoimmunity start...")
    list_of_proteins=autoimmunity(list_of_proteins, args.proteome1, args.working_dir, args.NERVE_dir, args.e_value, args.minlength, 
                                  args.mismatch, args.substitution)
    logging.debug("Done.")
    
    # Mouse immunity
    if mouse:
        logging.debug("Mouse immunity start...")
        list_of_proteins=mouse(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, args.proteome1,
                               args.minlength, args.substitution, args.mismatch)
        logging.debug("Done.")
        
    # Conservation
    #logging.debug(f'list of proteins before conservation:\n{list_of_proteins}')
    if args.proteome2:
        logging.debug("Conservation start...")
        list_of_proteins=conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, 
                                      args.proteome1, args.proteome2,
                                      args.minlength, args.substitution, args.mismatch)
        logging.debug("Done.")
        
    # Virulence
    #logging.debug(f'list of proteins before virulence:\n{list_of_proteins}')
    if args.virulent:
        logging.debug("Virulence start...")
        list_of_proteins=virulence(list_of_proteins, args.working_dir, args.iFeature_dir, args.proteome1, args.NERVE_dir)
        logging.debug("Done.")
        
    # annotation
    #logging.debug(f'list of proteins before annotation:\n{list_of_proteins}')
    if args.annotation:
        logging.debug("Annotation start...")
        list_of_proteins=annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        logging.debug("Done.")
    
    # select
    if args.select:
        logging.debug("Select start...")
        final_proteins=select(list_of_proteins, args.p_ad_no_citoplasm_filter, args.p_ad_extracellular_filter, 
               args.transmemb_doms_limit, args.padlimit, args.mouse, 
               args.mouse_peptides_sum_limit, args.virlimit, args.virulent)
        logging.debug("Done.")

    # final
    if args.proteome2:
        final_proteins.sort(key=lambda p: p.conservation_score, reverse=True) # ranking
    # return .csv outputs
    output(final_proteins, os.path.join(args.working_dir, 'vaccine_candidates.csv'))
    output(list_of_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'))
    #Protein.Protein.information_to_csv(list_of_proteins)
    logging.debug("Done: NERVE has finished its analysis!")
    
def bashCmdMethod(bashCmd):
    """Run bash commands"""
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

def quality_control(path_to_fasta:str, working_dir)->None:
    """
    Remove sequences with non-canonical aminoacid symbols. U (Se-Cys) is substituted with C (Cys). Returns
    "non_filtered_"+input_filename and overwrites input.
    param: path_to_fasta: full path to fasta file containing the proteome with .fasta extension;
    param: working_dir: working directory
    """
    def is_fasta(filename:str):
        """Function that rise an error if the format is not .fasta"""
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta) == True:
                return list(fasta)
            else:
                raise ValueError(f'{filename} is not in fasta format')
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    aa_dic = {'C': 'C', 'D': 'D', 'S': 'S', 'Q': 'Q', 'K': 'K', 'I': 'I', 'P': 'P', 'T': 'T', 'F': 'F', 'N': 'N', 
              'G': 'G', 'H': 'H', 'L': 'L', 'R': 'R', 'W': 'W', 'A': 'A', 'V': 'V', 'E': 'E', 'Y': 'Y', 'M': 'M', 
              'U':'C'}
    outlist = []
    #path_to_fasta = ".".join([path_to_fasta, ".fasta"])
    
    fasta_list = is_fasta(path_to_fasta)    
    filename = open(os.path.join(working_dir, "original_"+path_to_fasta.split('/')[-1]), 'w')
    SeqIO.write(fasta_list, filename, "fasta")
    filename.close()
    
    for record in fasta_list:
        flag = True
        new_seq =''
        for aa in str(record.seq):
            if aa not in aa_dic:
                flag = False
                logging.debug(f'Found non-canonical aminoacid named {aa} in sequence {record.id}')
            else:
                new_seq += aa_dic[aa]
        record.seq = Seq(new_seq)
        if flag == True:
            outlist.append(record)
        else:    
            logging.debug(f'Sequence {record.id} has been discarded for the presence of non-canonical aminoacids.')     
    # filtered
    filename = open(path_to_fasta, 'w')
    SeqIO.write(outlist, filename, "fasta")
    filename.close()
    return None

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
            if p.id in row["name"] and row["prediction"][0].reliability >= treshold:
                p.localization = row["prediction"][0].localization
    # save cello raw predictions
    df.to_csv(os.path.join(working_dir, 'cello_predictions.csv'))
    return list_of_proteins

def adhesin(list_of_proteins, working_dir, NERVE_dir)->list:
    "Runs adhesin predictions"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    # take the model from nerve but the methods from spaan, cause the model in spaan could be modified and tested
    model = keras.models.load_model(os.path.join(NERVE_dir, 'espaan_model.h5')) 
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
        

def razor(list_of_proteins, working_dir, transmemb_doms_limit, razor_len)->list:
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
            if len(longest_loop) > razor_len:
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
    
    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir, "sapiens_database/sapiens"), evalue=e_value, outfmt=5, out=os.path.join(working_dir,"sapiens.xml")) # 5 is for xml 
    stdout, stderr = blastx_cline()
    #logging.debug("Warning: you can find a sapiens.xml file on your working directory which is the outputs of the autoimmunity module.\nDo not delete during the computation!\nAfter the computation it will be deleted in order to avoid future collisions.")
    # for each result in the .xml file...
    outfile=open(os.path.join(working_dir, 'autoimmunity_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir,"sapiens.xml"))):
        query_name = record.query.split(' ')[0] # take only the query id 
        # take the right candidate to update
        tmp_protein = list_of_proteins[0]
        for p in list_of_proteins:
            if p.id == query_name:
                tmp_protein = p
        # for each effective alignment between the tmp candidate and the human proteome
        for alignment in record.alignments:
            for hsp in alignment.hsps: # collect all the interesting peptides
                #print(hsp.query, hsp.query_start, hsp.match)
                tmp_protein.list_of_shared_human_peps += Protein.Protein.hsp_match_parser(hsp.match, hsp.query, parsing_window_size=minlength, max_sub=substitution, max_mismatch=mismatch)
        # print out the peptides (if there are any)
        if len(tmp_protein.list_of_shared_human_peps) == 0:
            outfile.write("\nNo immunogenic peptides for " + query_name)   
        else:
            outfile.write("\nList of immunogenic peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    outfile.close()
    os.remove(os.path.join(working_dir, "sapiens.xml")) # delete after the computation
    
    # sum peptides
    logging.debug('Run sum of peptides')
    for p in list_of_proteins:
        p.sapiens_peptides_sum=0
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
    mhcpep = pandas.read_csv(os.path.join(NERVE_dir, "mhcpep/mhcpep_sapiens.csv"), skipinitialspace=True)
    number_of_proteins = len(list_of_proteins)
    for p in list_of_proteins:
        for seq in p.list_of_shared_human_peps:
            for pep in mhcpep['Epitope.2']:
                tmp_matches = Protein.Protein.peptide_comparison(seq, pep)
                p.list_of_peptides_from_comparison_with_mhcpep_sapiens += tmp_matches
    return list_of_proteins
            
def mouse(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, minlength, substitution, mismatch)->list:
    """Module to run mouse immunity check"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir,"mouse_database/mouse"), evalue=e_value, outfmt=5, out=os.path.join(working_dir, "mouse.xml"))
    stdout, stderr = blastx_cline()
    outfile=open(os.path.join(working_dir, 'mouse_immunity_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir, "mouse.xml"))):
        query_name = record.query.split(' ')[0]
        tmp_protein = list_of_proteins[0]
        for p in list_of_proteins:
            if p.id == query_name:
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
    
    mhcpep = pandas.read_csv(os.path.join(NERVE_dir, "mhcpep/mhcpep_mouse.csv"), skipinitialspace=True)
    number_of_proteins = len(list_of_proteins)
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
    
    for record in NCBIXML.parse(open(os.path.join(working_dir,"comparison.xml"))):
        query_name = record.query.split(' ')[0] 
    
        tmp_protein = list_of_proteins[0]
        for p in list_of_proteins:
            if p.id == query_name:
                tmp_protein = p
                #max_score = 0
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #if hsp.score > max_score: max_score = hsp.score
                tmp_protein.list_of_shared_conserv_proteome_peps += Protein.Protein.hsp_match_parser(hsp.match, hsp.query, parsing_window_size=minlength, max_sub=substitution, max_mismatch=mismatch)
    # sum peptides
    logging.debug('Run sum of peptides')
    for p in list_of_proteins:
            p.conservation_score = 0
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
    virulent_model = tensorflow.keras.models.load_model(os.path.join(NERVE_dir, 'virulent_classification_model.h5'))
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
            if p.id in row['Protein']:
                p.annotations = row['Function']
    return list_of_proteins

def select(list_of_proteins, p_ad_no_citoplasm_filter, p_ad_extracellular_filter, transmemb_doms_limit,
           padlimit, mouse, mouse_peptides_sum_limit, virlimit, virulent)->list:
    """Selects suitable candidate proteins for vaccine production"""
    final_list = []
    for protein in list_of_proteins:
        if protein.localization == "Cytoplasmic": continue 
        if protein.p_ad < p_ad_no_citoplasm_filter and not protein.localization == "Extracellular": continue 
        if protein.p_ad < p_ad_extracellular_filter and protein.localization == "Extracellular": continue 
        if (protein.transmembrane_doms >= transmemb_doms_limit) and (protein.original_sequence_if_razor is None): continue
        if protein.sapiens_peptides_sum > .15: continue
        if len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens) >= 1: continue
        if (protein.localization == "Unknown") and (protein.p_ad < padlimit): continue
        if mouse:
            if protein.mouse_peptides_sum > mouse_peptides_sum_limit: continue 
            if len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse) >= 1: continue 
        if virulent:
            if protein.p_vir < virlimit: continue
        final_list.append(protein)
    return final_list

def output(list_of_proteins, outfile):
    pd.DataFrame([[str(protein.id),
                 str(protein.accession),
                 str(protein.sequence),
                 str(protein.original_sequence_if_razor),
                 str(protein.length),
                 str(protein.localization),
                 str(protein.p_ad),
                 str(protein.transmembrane_doms),
                 str(protein.tmhmm_seq),
                 str(protein.list_of_shared_human_peps),
                 str(protein.list_of_shared_mouse_peps),
                 str(protein.list_of_shared_conserv_proteome_peps),
                 str(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens),
                 str(protein.list_of_peptides_from_comparison_with_mhcpep_mouse),
                 str(protein.razor_loops),
                 str(protein.p_vir),
                 str(protein.sapiens_peptides_sum),
                 str(protein.mouse_peptides_sum),
                 str(protein.conservation_score),
                 str(protein.annotations)
                 ] for protein in list_of_proteins
                ], 
                columns= ['id ',
                    'accession',
                    'sequence',
                    'original_sequence_if_razor',
                    'length',
                    'localization',
                    'p_ad',
                    'transmembrane_doms',
                    'tmhmm_seq',
                    'list_of_shared_human_peps',
                    'list_of_shared_mouse_peps',
                    'list_of_shared_conserv_proteome_peps',
                    'list_of_peptides_from_comparison_with_mhcpep_sapiens',
                    'list_of_peptides_from_comparison_with_mhcpep_mouse',
                    'razor_loops',
                    'p_vir',
                    'sapiens_peptides_sum',
                    'mouse_peptides_sum',
                    'conservation_score',
                    'annotations'
                     ]
                ).to_csv(outfile) 
    
if __name__ == "__main__":
    main()
