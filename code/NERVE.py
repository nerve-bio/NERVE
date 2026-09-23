#!/usr/local/bin/python
"""Run NERVE, reverse vaccinology software"""

import argparse, logging, math, os, pickle, sys, time
from typing import NamedTuple

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # disable warnings

# Keep numerical predictions reproducible.
for _thread_var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS',
                    'TF_NUM_INTRAOP_THREADS', 'TF_NUM_INTEROP_THREADS'):
    os.environ[_thread_var] = '1'

try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    pass

from code.Protein import Protein
from code.Utils import dir_path
from code.Function import annotation
from code.Adhesin import extract_features, adhesin_predict
from code.Virulent_factor import virulent_factor_predict
from code.Quality_control import proteome_downloader, quality_control
from code.Subcellular import psortb
from code.Topology import tmhelices
from code.Razor import razor
from code.Immunity import  autoimmunity, conservation, mouse
from code.Select import output, select, scorer
from code.Epitope import epitope
from code.Checkpoint import Checkpointer, compute_signature

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
    ep_plots:bool
    threads:int
    keep_snapshots:str
    epitope_workers:int
    checkpoint:str
    resume:str
    max_unknown_fraction:float
    
    def print_args(self):
        return (f'''annotation: {self.annotation}, e_value: {self.e_value}, gram: {self.gram},
                minlength: {self.minlength}, mismatch: {self.mismatch}, mouse: {self.mouse},
                mouse_peptides_sum_limit: {self.mouse_peptides_sum_limit}, proteome1: {self.proteome1},
                proteome2: {self.proteome2}, padlimit: {self.padlimit}, razor: {self.razor}, razlen: {self.razlen}, select: {self.select},
                substitution: {self.substitution}, transmemb_doms_limit: {self.transmemb_doms_limit},
                virlimit: {self.virlimit}, virulent: {self.virulent},epitopes: {self.epitopes}, mhci_length: {self.mhci_length},
                mhcii_length: {self.mhcii_length}, mhci_overlap: {self.mhci_overlap}, mhcii_overlap: {self.mhcii_overlap},
                epitope_percentile: {self.epitope_percentile}, working_dir: {self.working_dir},
                NERVE_dir: {self.NERVE_dir}, iFeature_dir: {self.iFeature_dir},  DeepFri_dir: {self.DeepFri_dir}, ep_plots: {self.ep_plots}''')
    
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
                        default="False",
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
    parser.add_argument('-pl','--padlimit',
                        metavar='\b',
                        help="Set the probability of adhesin (pad) value cut-off for all proteins in the select module. Thus, proteins with a pad value < cut-off are discarded (0.-1)",
                        type=float,
                        default=0.5,
                        required=False,
                        )
    parser.add_argument('-rz','--razor',
                        metavar='\b', 
                        help="Activation (True) or deactivation (False) of the loop-razor module. This module allows the recovery of protein vaccine candidates, with more than 2 transmembrane domains, that would otherwise be discarded in the select module. The longest loop with minimum len == 'razlen' aa will replace the original protein sequence for following NERVE steps",
                        type=str,
                        default="False",
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
                        help="Parameter of select module. Proteins with trasmembrane domains >= transmemb_doms_limit are discarded",
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
                        default="False",
                        required=False,
                        )
    parser.add_argument('-ep', '--epitopes',
                        metavar='\b',
                        type=str,
                        help='Activate or deactivate epitope prediction module',
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
                        help='percentile decision threshold on which to predict epitopes from full length proteins',
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
    parser.add_argument('-epp', '--ep_plots',
                        metavar='\b',
                        help='Epitopes plotting script',
                        type=str,
                        required=False,
                        default="True"
                       )
    parser.add_argument('-t', '--threads',
                        metavar='\b',
                        help='Number of threads used by blastp in the immunity modules',
                        type=int,
                        required=False,
                        default=1
                       )
    parser.add_argument('-ks', '--keep_snapshots',
                        metavar='\b',
                        help="Pickle the protein list after every stage into 'snapshots/' for "
                             "comparing two runs stage by stage (see harness/)",
                        type=str,
                        required=False,
                        default="False"
                       )
    parser.add_argument('-ew', '--epitope_workers', type=int, default=0,
                        help='Epitope worker processes; 0 auto-detects CPU/cgroup capacity')
    parser.add_argument('-ck', '--checkpoint', type=str, default='False',
                        help='Write an atomic restart checkpoint after every completed stage')
    parser.add_argument('-rs', '--resume', type=str, default='False',
                        help='Resume from a compatible checkpoint (also keeps checkpointing enabled)')
    parser.add_argument('--max_unknown_fraction', type=float, default=0.05,
                        help='QC discard threshold for residues still unknown after U/O/B/Z/J mapping')


    args = parser.parse_args()

    for name in ('annotation', 'mouse', 'razor', 'select', 'virulent', 'epitopes',
                 'ep_plots', 'keep_snapshots', 'checkpoint', 'resume'):
        if getattr(args, name) not in ('True', 'False'):
            parser.error(f'--{name} must be True or False')
    if args.gram not in ('n', 'p'):
        parser.error('--gram must be n or p')
    for name in ('threads', 'minlength', 'razlen', 'transmemb_doms_limit'):
        if getattr(args, name) < 1:
            parser.error(f'--{name} must be positive')
    if args.epitope_workers < 0:
        parser.error('--epitope_workers must be 0 (automatic) or positive')
    for name in ('mismatch', 'substitution'):
        if not 0 <= getattr(args, name) <= args.minlength:
            parser.error(f'--{name} must be between 0 and --minlength')
    for name in ('padlimit', 'virlimit', 'epitope_percentile',
                 'max_unknown_fraction'):
        value = getattr(args, name)
        if not math.isfinite(value) or not 0 <= value <= 1:
            parser.error(f'--{name} must be a finite number between 0 and 1')
    if not math.isfinite(args.mouse_peptides_sum_limit) or args.mouse_peptides_sum_limit < 0:
        parser.error('--mouse_peptides_sum_limit must be a finite nonnegative number')
    if args.mouse == 'True' and args.mouse_peptides_sum_limit == 0:
        parser.error('--mouse_peptides_sum_limit must be positive when --mouse True')
    if not math.isfinite(args.e_value) or args.e_value <= 0:
        parser.error('--e_value must be a finite positive number')

    return Args(args.annotation, args.e_value, args.gram, args.minlength, args.mismatch,
                args.mouse, args.mouse_peptides_sum_limit, args.proteome1, args.proteome2,
                args.padlimit, args.razor,
                args.razlen, args.select, args.substitution, args.transmemb_doms_limit, args.virlimit,
                args.virulent, args.epitopes,
                args.mhci_length, args.mhcii_length, args.mhci_overlap, args.mhcii_overlap,
                args.epitope_percentile,args.working_dir, args.NERVE_dir, args.iFeature_dir, args.DeepFri_dir, args.ep_plots,
                args.threads, args.keep_snapshots, args.epitope_workers, args.checkpoint,
                args.resume, args.max_unknown_fraction)


# Per-stage timings.
_STAGE_TIMINGS = []
_CHECKPOINTER = None


def record_stage(args, ordinal, name, proteins, started, all_proteins=None):
    """Record a completed pipeline stage."""
    if _CHECKPOINTER is not None and _CHECKPOINTER.done(ordinal):
        return
    _STAGE_TIMINGS.append((ordinal, name, time.monotonic() - started))
    if _CHECKPOINTER is not None:
        _CHECKPOINTER.save(ordinal, name, {
            'args': args,
            'proteins': proteins,
            'all_proteins': all_proteins if all_proteins is not None else proteins,
        })
    if str(getattr(args, 'keep_snapshots', 'False')) != "True":
        return
    directory = os.path.join(args.working_dir, 'snapshots')
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f'{ordinal:02d}_{name}.pkl')
    # Replace the snapshot atomically.
    with open(path + '.tmp', 'wb') as handle:
        pickle.dump(proteins, handle)
    os.replace(path + '.tmp', path)


def timing_report() -> str:
    """Format the per-stage timings as a table with totals and percentages."""
    if not _STAGE_TIMINGS:
        return ''
    total = sum(seconds for _, _, seconds in _STAGE_TIMINGS)
    lines = ['', '=' * 50, '{:^50}'.format('Per-stage timing'), '=' * 50,
             f'{"#":>2}  {"stage":<20}{"seconds":>10}{"%":>8}']
    for ordinal, name, seconds in sorted(_STAGE_TIMINGS):
        share = (100 * seconds / total) if total else 0
        lines.append(f'{ordinal:>2}  {name:<20}{seconds:>10.1f}{share:>7.1f}%')
    lines.append('-' * 40)
    lines.append(f'{"":>2}  {"TOTAL":<20}{total:>10.1f}{100.0:>7.1f}%')
    return '\n'.join(lines)


def main():
    """Runs NERVE"""
    # record time:
    nerve_start = time.monotonic()
    global _CHECKPOINTER
    args=get_args()
    print("Start NERVE 2.0")
    
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
            raise RuntimeError(f'Unable to download proteome1 {args.proteome1}: {e}') from e
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
            except Exception as e:
                raise RuntimeError(f'Unable to download proteome2 {args.proteome2}: {e}') from e
            logging.debug(f'{args.proteome2} successfully downloaded')
            args = args._replace(proteome2=os.path.join(args.working_dir,'proteome2.fasta'))   
    signature_fields = args._asdict().copy()
    for field in ('working_dir', 'proteome1', 'proteome2', 'keep_snapshots',
                  'checkpoint', 'resume'):
        signature_fields.pop(field, None)
    source_identity = {
        'commit': os.environ.get('NERVE_SOURCE_COMMIT', 'unknown'),
        'tree': os.environ.get('NERVE_SOURCE_TREE', 'unknown'),
        'dirty': os.environ.get('NERVE_SOURCE_DIRTY', 'unknown'),
    }
    signature = compute_signature(
        [args.proteome1, args.proteome2], signature_fields, source_identity)
    _CHECKPOINTER = Checkpointer(
        args.working_dir, signature,
        enabled=args.checkpoint == 'True', resume=args.resume == 'True')
    restored = _CHECKPOINTER.restored_state()
    completed_stage = _CHECKPOINTER.bundle['ordinal'] if restored else 0
    if restored:
        args = restored['args']
        list_of_proteins = restored['all_proteins']
        final_proteins = restored['proteins']
        print(f'Resuming after stage {completed_stage}')
    logging.info('NERVE source commit=%s tree=%s dirty=%s run_signature=%s',
                 source_identity['commit'], source_identity['tree'],
                 source_identity['dirty'], signature)
    print("10% done")
    
    # run quality control
    if completed_stage < 1:
        start=time.monotonic()
        logging.debug(f'Start quality control of proteome1 ({args.proteome1})')
        list_of_fasta_proteins, proteome1_new_path = quality_control(
            args.proteome1, args.working_dir, upload=True,
            max_unknown_fraction=args.max_unknown_fraction)
        args=args._replace(proteome1=proteome1_new_path)
        if len(list_of_fasta_proteins) == 0:
            raise ValueError('All input protein sequences have been discarded.')
        if args.proteome2:
            proteome2_new_path=quality_control(
                args.proteome2, args.working_dir,
                max_unknown_fraction=args.max_unknown_fraction)
            args = args._replace(proteome2=proteome2_new_path)
        list_of_proteins = []
        for p in list_of_fasta_proteins:
            protein = Protein(str(p.description), str(p.seq))
            protein.original_id = str(getattr(p, 'original_description', p.description))
            protein.accession = str(getattr(p, 'accession', protein.accession))
            list_of_proteins.append(protein)
        record_stage(args, 1, 'quality_control', list_of_proteins, start)
    
    # subcellular localization prediction
    start = time.monotonic()
    logging.debug("Subcelloc start with psortb...")
    if completed_stage < 2:
        list_of_proteins = psortb(list_of_proteins, args.working_dir, args.gram, args.proteome1)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end - start))
    print("20% done")
    record_stage(args, 2, 'localization', list_of_proteins, start)

    # Adhesin
    logging.debug("Adhesin start...")
    start=time.monotonic()
    # extract features
    if completed_stage < 3:
        list_of_proteins = extract_features(list_of_proteins, args.NERVE_dir, args.iFeature_dir, args.working_dir, args.proteome1)
        list_of_proteins = adhesin_predict(list_of_proteins, args.NERVE_dir)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("30% done")
    record_stage(args, 3, 'adhesin', list_of_proteins, start)
    
    # Tmhelices
    logging.debug("Tmhelices start...")
    start=time.monotonic()
    if completed_stage < 4:
        list_of_proteins = tmhelices(list_of_proteins, args.working_dir,
                                     threads=args.threads)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("40% done")
    record_stage(args, 4, 'topology', list_of_proteins, start)
    
    # Razor
    if args.razor == "True":
        logging.debug("Loop-razor start...")
        start=time.monotonic()
        if completed_stage < 5:
            list_of_proteins=razor(list_of_proteins, args.working_dir, args.transmemb_doms_limit, args.razlen)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
        record_stage(args, 5, 'razor', list_of_proteins, start)
    print("50% done")
    
    # Autoimmunity
    logging.debug("Autoimmunity start...")
    start=time.monotonic()
    if completed_stage < 6:
        list_of_proteins=autoimmunity(list_of_proteins, args.proteome1, args.working_dir, args.NERVE_dir, args.e_value,
                                      args.minlength, args.mismatch, args.substitution, threads=args.threads)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    print("60% done")
    record_stage(args, 6, 'autoimmunity', list_of_proteins, start)
    
    # Mouse immunity
    if args.mouse == "True":
        start = time.monotonic()
        logging.debug("Mouse immunity start...")
        if completed_stage < 7:
            list_of_proteins = mouse(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, args.proteome1,
                                     args.minlength, args.substitution, args.mismatch, threads=args.threads)
        end = time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 7, 'mouse', list_of_proteins, start)
    
    # Conservation
    #logging.debug(f'list of proteins before conservation:\n{list_of_proteins}')
    if args.proteome2:
        start = time.monotonic()
        logging.debug("Conservation start...")
        if completed_stage < 8:
            list_of_proteins = conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value,
                                            args.proteome1, args.proteome2,
                                            args.minlength, args.substitution, args.mismatch, threads=args.threads)
        end = time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 8, 'conservation', list_of_proteins, start)
    print("70% done")
        
    # Virulence
    #logging.debug(f'list of proteins before virulence:\n{list_of_proteins}')
    if args.virulent == "True":
        start=time.monotonic()
        logging.debug("Virulence start...")
        if completed_stage < 9:
            list_of_proteins = virulent_factor_predict(list_of_proteins, args.NERVE_dir)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 9, 'virulence', list_of_proteins, start)
    print("80% done")
        
    # annotation
    if args.annotation == "True":
        start = time.monotonic()
        logging.debug("Annotation start...")
        if completed_stage < 10:
            list_of_proteins = annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 10, 'annotation', list_of_proteins, start)
    print("90% done")
    
    # score
    for protein in list_of_proteins:
        protein.score = scorer(protein, args.mouse_peptides_sum_limit, args.mouse)
    
    # select
    if completed_stage < 11:
        final_proteins=list_of_proteins
    if args.select == "True":
        logging.debug("Select start...")
        start=time.monotonic()
        if completed_stage < 11:
            final_proteins = select(list_of_proteins, args.transmemb_doms_limit,
                                    args.padlimit, args.mouse, args.mouse_peptides_sum_limit, args.virlimit, args.virulent, args.razor)
        end = time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 11, 'select', final_proteins, start, list_of_proteins)

    #if args.virulent == "True":
    #    final_proteins.sort(key = lambda p: p.p_vir, reverse = True)
    
    # 12.Epitope prediction
    if args.epitopes == "True":
        #print("=" * 50)
        #print("{:^50}".format('Epitope prediction of best candidates with epitopepredict starts'))
        #print("=" * 50)
        start = time.monotonic()
        logging.debug('Epitope prediction starts ...')
        if completed_stage < 12:
            final_proteins = epitope(final_proteins,
                                     args.working_dir, args.mhci_length, args.mhcii_length,
                                     args.mhci_overlap, args.mhcii_overlap, args.epitope_percentile,
                                     args.ep_plots, args.transmemb_doms_limit, args.epitope_workers)
        end = time.monotonic()
        logging.debug(f'Epitope prediction done in {end - start} seconds')
        record_stage(args, 12, 'epitope', final_proteins, start, list_of_proteins)
    
    # return .csv outputs
    output(final_proteins, os.path.join(args.working_dir, 'vaccine_candidates.csv'))
    # collect discarded proteins
    final_proteins_names = [p.id for p in final_proteins]
    discarded_proteins = [p for p in list_of_proteins if p.id not in final_proteins_names]
    output(discarded_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'))
    
    nerve_end = time.monotonic()
    report = timing_report()
    if report:
        print(report)
        logging.debug(report)
    logging.debug("Done: NERVE has finished its analysis in: {:.4f} seconds".format(nerve_end-nerve_start))
    print("100% done")
    print("End NERVE computation successfully.")
    _CHECKPOINTER.finish()

if __name__ == "__main__":
    main()
