#!/usr/local/bin/python
"""Run NERVE, reverse vaccinology software"""

import argparse, logging, os, time
from typing import NamedTuple

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # disable warnings

from Protein import *
from Utils import bashCmdMethod, dir_path
from Function import annotation
from Adhesin import extract_features, adhesin_predict
from Virulent_factor import virulent_factor_predict
from Quality_control import proteome_downloader, proteome_uploader, quality_control
from Subcellular import psortb
from Topology import tmhelices
from Razor import razor
from Immunity import  autoimmunity, conservation, mouse
from Select import output, select
from Epitope import *

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
    transmemb_doms_limit:int
    virlimit:float
    virulent:bool
    epitopes:str
    mhci_length:int
    mhcii_length:int
    mhci_overlap:int
    mhcii_overlap:int
    epitope_percentile:float
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
                substitution: {self.substitution}, transmemb_doms_limit: {self.transmemb_doms_limit},
                virlimit: {self.virlimit}, virulent: {self.virulent},epitopes: {self.epitopes}, mhci_length: {self.mhci_length},
                mhcii_length: {self.mhcii_length}, mhci_overlap: {self.mhci_overlap}, mhcii_overlap: {self.mhcii_overlap},
                epitope_percentile: {self.epitope_percentile}, working_dir: {self.working_dir},
                NERVE_dir: {self.NERVE_dir}, iFeature_dir: {self.iFeature_dir},  DeepFri_dir: {self.DeepFri_dir}''')
    
def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description="Run vaccine candidate prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-a','--annotation',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of annotation module. Uses DeepFri to retrieve protein functional onthologies",
                        type=str,
                        required=False,
                        default="True"
                        )
    parser.add_argument('-ev','--e_value',
                        metavar='\b', 
                        help="Expect-value used in blastp for immunity modules",
                        type=float,
                        default=1e-10,
                        required=False,
                        )
    parser.add_argument('-g','--gram',
                        metavar='\b', 
                        help="Negative (n) or positive (p) gram stain of the pathogen of interest",
                        type=str,
                        required=True,
                        )
    parser.add_argument('-ml','--minlength',
                        metavar='\b', 
                        help="Minimal length required for shared peptides to be extracted in comparison analyses versus human and/or mouse",
                        type=int,
                        default=9,
                        required=False,
                        )
    parser.add_argument('-mm','--mismatch',
                        metavar='\b', 
                        help="Maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules",
                        type=int,
                        default=1,
                        required=False,
                        )
    parser.add_argument('-m','--mouse',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the eventual shared peptides is carried out as in the autoimmunity module",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-mpsl','--mouse_peptides_sum_limit',
                        metavar='\b', 
                        help="Parameter calculated in mouse module and used by select module. Protein with 'sum of shared peptides of the i-protein with mouse proteins/number of aminoacids of the i-protein' <= mouse_peptides_sum_limit and with absence of match mhc-I and Mhc-II mouse ligands are selected",
                        type=float,
                        default=0.15,
                        required=False,
                        )
    parser.add_argument('-p1','--proteome1',
                        metavar='\b', 
                        help='Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score)',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-p2','--proteome2',
                        metavar='\b', 
                        help='Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score)',
                        type=str,
                        required=False,
                        )
    parser.add_argument('-paefilter','--p_ad_extracellular_filter',
                        metavar='\b', 
                        help="Parameter of select module. Extracellular proteins with a probability of adhesin (pad) lower than p_ad_extracellular_filter are discarded (0.-1)",
                        type=float,
                        default=0.38,
                        required=False,
                        )
    parser.add_argument('-pacfilter','--p_ad_no_citoplasm_filter',
                        metavar='\b', 
                        help="Parameter of select module. Non-cytoplasmic Proteins with a probability of adhesin (pad) lower than p_ad_no_citoplasm_filter are discarded (0.-1)",
                        type=float,
                        default=0.46,
                        required=False,
                        )    
    parser.add_argument('-pl','--padlimit',
                        metavar='\b',
                        help="Set the probability of adhesin (pad) value cut-off for proteins with 'Unknown' localization in the select module. Thus, these proteins with a pad value < cut-off are discarded (0.-1)",
                        type=float,
                        default=0.85,
                        required=False,
                        )
    parser.add_argument('-rz','--razor',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of the loop-razor module. This module allows the recovery of protein vaccine candidates, with more than 2 transmembrane domains, that would otherwise be discarded in the select module. The longest loop with minimum len == 'razlen' aa will replace the original protein sequence for following NERVE steps",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-rl','--razlen',
                        metavar='\b', 
                        help="Set minimal length of loop considered in loop-razor module",
                        type=int,
                        default=50,
                        required=False,
                        )
    parser.add_argument('-s','--select',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of select module, which filters PVC from proteome1",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-ss','--substitution',
                        metavar='\b', 
                        help="Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-tdl','--transmemb_doms_limit',
                        metavar='\b', 
                        help="Parameter of select module. Proiteins with trasmembrane domains >= transmemb_doms_limit are discarded",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-vl','--virlimit',
                        metavar='\b', 
                        help="Cut-off value for NERVirulent in the select module (0.-1)",
                        type=float,
                        default=0.5,
                        required=False,
                        )    
    parser.add_argument('-vir','--virulent',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of NERVirulent module, predictor of the probability of being a virulence factor",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-ep', '--epitopes',
                        metavar='\b',
                        type=str,
                        help='Activate or deactivate epitopes module',
                        required=False,
                        default="True"
                        )
    parser.add_argument('-m1l', '--mhci_length',
                        metavar='\b',
                        type=int,
                        help='mhci binders length (9, 10, 11 are available)',
                        required=False,
                        choices=[9, 10, 11],
                        default=9
                        )
    parser.add_argument('-m2l', '--mhcii_length',
                        metavar='\b',
                        type=int,
                        help='mhcii binders length (9, 11, 13, 15 are available)',
                        required=False,
                        choices=[9, 11, 13, 15],
                        default=11
                        )
    parser.add_argument('-m1ovr', '--mhci_overlap',
                        metavar='\b',
                        type=int,
                        help='mhci-epitope overlap',
                        required=False,
                        choices=[1, 2],
                        default=1
                        )
    parser.add_argument('-m2ovr', '--mhcii_overlap',
                        metavar='\b',
                        type=int,
                        help='mhcii-epitope overlap',
                        required=False,
                        choices=[1, 2],
                        default=1
                        )
    parser.add_argument('-prt', '--epitope_percentile',
                        metavar='\b',
                        type=float,
                        help='percentile decision threshold on whick to predict epitopes from full length proteins',
                        required=False,
                        default=0.9
                        )
     
    parser.add_argument('-wd','--working_dir',
                        metavar='\b', 
                        help='Path to working directory. If not existing, a working directory with the given path is created',
                        type=str,
                        required=False,
                        default='./'
                        )
    parser.add_argument('-nd','--NERVE_dir',
                        metavar='\b', 
                        help='Path to NERVE repository folder (download from: https://github.com/nicolagulmini/NERVE)',
                        type=dir_path,
                        required=False,
                        default='/usr/nerve_python/NERVE'
                        )
    parser.add_argument('-id','--iFeature_dir',
                        metavar='\b', 
                        help='Path to iFeature repository folder (download from: https://github.com/Superzchen/iFeature)',
                        type=dir_path,
                        required=False,
                        default='/usr/nerve_python/assets/iFeature'
                        )
    parser.add_argument('-dfd','--DeepFri_dir',
                        metavar='\b', 
                        help='Path to DeepFri folder (download from: https://github.com/flatironinstitute/DeepFRI)',
                        type=dir_path,
                        required=False,
                        default='/usr/nerve_python/assets/DeepFri'
                        )
    
    
    args = parser.parse_args()

    return Args(args.annotation, args.e_value, args.gram, args.minlength, args.mismatch,
                args.mouse, args.mouse_peptides_sum_limit, args.proteome1, args.proteome2, 
                args.p_ad_extracellular_filter, args.p_ad_no_citoplasm_filter, args.padlimit, args.razor, 
                args.razlen, args.select, args.substitution, args.transmemb_doms_limit, args.virlimit, 
                args.virulent, args.epitopes,
                args.mhci_length, args.mhcii_length, args.mhci_overlap, args.mhcii_overlap,
                args.epitope_percentile,args.working_dir, args.NERVE_dir, args.iFeature_dir, args.DeepFri_dir)


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
    if os.path.isdir(args.working_dir) == False:
        os.makedirs(args.working_dir)
    # define log file
    logging.basicConfig(filename=os.path.join(args.working_dir, 'logfile.log'),
                        filemode='w',
                        level=logging.DEBUG,
                        force=True)
    logging.debug(f'Running NERVE with the following parameters:\n{args.print_args()}')    
    # check input and download proteome:
    if os.path.isfile(args.proteome1) == True:
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
        if os.path.isfile(args.proteome2) == True:
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
    list_of_fasta_proteins, proteome1_new_path = quality_control(args.proteome1, args.working_dir, upload=True)
    # update input path of proteome1
    args=args._replace(proteome1=proteome1_new_path)
    logging.debug(f'Finish quality control of proteome1. Updated path: ({args.proteome1})')
    if len(list_of_fasta_proteins) == 0:
        raise ValueError(f'All input protein sequences have been discarded. See {os.path.join(args.working_dir, "logfile.log")} for more information.')
    if args.proteome2:
        logging.debug(f'Start quality control of proteome2 ({args.proteome2})')
        proteome2_new_path=quality_control(args.proteome2, args.working_dir)
        # update proteome2 new path 
        args = args._replace(proteome2=proteome2_new_path)
        logging.debug(f'Finish quality control of proteome2. Updated path: ({args.proteome2})')
    logging.debug(f'Extract protein sequences and IDs from proteome1')
    list_of_proteins = []
    for p in list_of_fasta_proteins:
        p_id = str(p.name)
        p_seq = str(p.seq)
        list_of_proteins.append(Protein(p_id, p_seq))
    end=time.time()
    logging.debug(f'{len(list_of_fasta_proteins)} proteins loaded in {end-start} seconds')
            
    # subcellular localization prediction
    start = time.time()
    logging.debug("Subcelloc start with psortb...")
    list_of_proteins = psortb(list_of_proteins, args.working_dir, args.gram, args.proteome1)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end - start))
    print("20% done")
    
    # Adhesin
    logging.debug("Adhesin start...")
    start=time.time()
    # extract features
    list_of_proteins = extract_features(list_of_proteins, args.NERVE_dir, args.iFeature_dir, args.working_dir, args.proteome1)
    list_of_proteins = adhesin_predict(list_of_proteins, args.NERVE_dir)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("30% done")
    
    # Tmhelices
    logging.debug("Tmhelices start...")
    start=time.time()
    list_of_proteins = tmhelices(list_of_proteins, args.working_dir)
    end=time.time()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("40% done")
    
    # Razor
    if args.razor == "True":
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
    if args.mouse == "True":
        start = time.time()
        logging.debug("Mouse immunity start...")
        list_of_proteins = mouse(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, args.proteome1,
                               args.minlength, args.substitution, args.mismatch)
        end = time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
    
    # Conservation
    #logging.debug(f'list of proteins before conservation:\n{list_of_proteins}')
    if args.proteome2:
        start = time.time()
        logging.debug("Conservation start...")
        list_of_proteins = conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, 
                                      args.proteome1, args.proteome2,
                                      args.minlength, args.substitution, args.mismatch)
        end = time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
    print("70% done")
        
    # Virulence
    #logging.debug(f'list of proteins before virulence:\n{list_of_proteins}')
    if args.virulent == "True":
        start=time.time()
        logging.debug("Virulence start...")
        list_of_proteins = virulent_factor_predict(list_of_proteins, args.NERVE_dir)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
    print("80% done")
        
    # annotation
    #logging.debug(f'list of proteins before annotation:\n{list_of_proteins}')
    if args.annotation == "True":
        start = time.time()
        logging.debug("Annotation start...")
        list_of_proteins = annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        end=time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
    print("90% done")
    
    # select
    final_proteins=list_of_proteins
    if args.select == "True":
        logging.debug("Select start...")
        start=time.time()
        final_proteins = select(list_of_proteins, args.transmemb_doms_limit,
                                args.padlimit, args.mouse, args.mouse_peptides_sum_limit, args.virlimit, args.virulent, args.annotation)
        end = time.time()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))

    #if args.virulent == "True":
    #    final_proteins.sort(key = lambda p: p.p_vir, reverse = True)
    
    # 12.Epitope prediction
    if args.epitopes == "True":
        print("=" * 50)
        print("{:^50}".format('Epitope prediction of best candidates with epitopepredict starts'))
        print("=" * 50)
        start = time.time()
        logging.debug('Epitope prediction starts ...')
        final_proteins = epitope(final_proteins, args.mouse, args.mouse_peptides_sum_limit,
                                 args.working_dir, args.mhci_length, args.mhcii_length,
                                 args.mhci_overlap, args.mhcii_overlap, args.epitope_percentile)
        end = time.time()
        logging.debug(f'Epitope prediction done in {end - start} seconds')
    
    # return .csv outputs
    output(final_proteins, os.path.join(args.working_dir, 'vaccine_candidates.csv'), args.mouse_peptides_sum_limit, args.mouse)
    # collect discarded proteins
    final_proteins_names = [p.id for p in final_proteins]
    discarded_proteins = [p for p in list_of_proteins if p.id not in final_proteins_names]
    output(discarded_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'), args.mouse_peptides_sum_limit, args.mouse)
    
    nerve_end = time.time()
    logging.debug("Done: NERVE has finished its analysis in: {:.4f} seconds".format(nerve_end-nerve_start))
    print("100% done")
    print("End NERVE computation successfully.")

if __name__ == "__main__":
    main()
