#!/usr/local/bin/python
"""Run NERVE, reverse vaccinology software"""

import argparse, logging, math, os, pickle, signal, sys, time
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
from code.Select import output, select, scorer, write_selection_reasons
from code.Epitope import epitope
from code.Checkpoint import Checkpointer, compute_signature
from code.Run_record import RunRecord, file_hash, input_records

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
    debug:str
    
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
                        help="Pickle the protein list after every stage into 'snapshots/'",
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
    parser.add_argument('--debug', choices=['True', 'False'], default='False',
                        help='Write detailed developer diagnostics to debug.log')


    args = parser.parse_args()

    for name in ('annotation', 'mouse', 'razor', 'select', 'virulent', 'epitopes',
                 'ep_plots', 'keep_snapshots', 'checkpoint', 'resume', 'debug'):
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
                args.resume, args.max_unknown_fraction, args.debug)


# Per-stage timings.
_STAGE_TIMINGS = []
_CHECKPOINTER = None
_RUN_RECORD = None
_RUN_STARTED = None


def check_output_directory(directory, resume):
    """Reject occupied output directories."""
    path = os.path.abspath(directory)
    checkpoint = os.path.join(path, 'checkpoint.pkl')
    manifest = os.path.join(path, 'run_manifest.json')
    if resume:
        if not os.path.isfile(checkpoint) or not os.path.isfile(manifest):
            raise ValueError(f'Cannot resume {path}: compatible checkpoint and manifest required')
        return
    if not os.path.isdir(path):
        return
    occupied = {'run_manifest.json', 'run.log', 'debug.log', 'logfile.log',
                'checkpoint.pkl', 'vaccine_candidates.csv', 'discarded_proteins.csv',
                'selection_reasons.csv', 'conservation_raw_output.txt',
                'comparison.xml', 'sapiens.xml', 'mouse.xml',
                '_MF_predictions.csv', '_MF_pred_scores.json',
                'compare_proteome',
                'epitope', 'snapshots'}
    found = [name for name in os.listdir(path)
             if name in occupied or name.startswith(('cleaned_', 'discarded_sequences_',
                                                     'qc_discard_log_', 'qc_modifications_log_'))]
    if found:
        raise ValueError(f'Refusing fresh run in {path}: existing result {found[0]}; '
                         'choose a new directory or use --resume True with a checkpoint')


def configure_logging(directory, debug, resume):
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)
        handler.close()
    root.setLevel(logging.DEBUG if debug else logging.INFO)
    formatter = logging.Formatter('%(message)s')
    for handler in (logging.StreamHandler(sys.stdout),
                    logging.FileHandler(os.path.join(directory, 'run.log'),
                                        mode='a' if resume else 'w')):
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        root.addHandler(handler)
    if debug:
        handler = logging.FileHandler(os.path.join(directory, 'debug.log'),
                                      mode='a' if resume else 'w')
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s'))
        root.addHandler(handler)


def start_stage(ordinal, name, count, completed_stage):
    if completed_stage < ordinal:
        logging.info('Stage %s %s: starting (%s proteins)', ordinal, name, count)


def qc_counts(input_path, directory):
    base = os.path.basename(input_path)
    stem = os.path.splitext(base)[0]
    cleaned = os.path.join(directory, 'cleaned_' + base)
    discarded = os.path.join(directory, f'qc_discard_log_{stem}.tsv')
    modified = os.path.join(directory, f'qc_modifications_log_{stem}.tsv')
    with open(cleaned) as handle:
        kept = sum(line.startswith('>') for line in handle)
    rejected = 0
    if os.path.isfile(discarded):
        with open(discarded) as handle:
            rejected = max(0, sum(1 for _ in handle) - 1)
    edited = set()
    if os.path.isfile(modified):
        with open(modified) as handle:
            next(handle, None)
            edited = {line.split('\t', 1)[0] for line in handle}
    return {'kept': kept, 'discarded': rejected, 'modified': len(edited)}


def record_stage(args, ordinal, name, proteins, started, all_proteins=None):
    """Record a completed pipeline stage."""
    if _CHECKPOINTER is not None and _CHECKPOINTER.done(ordinal):
        return
    seconds = time.monotonic() - started
    _STAGE_TIMINGS.append((ordinal, name, seconds))
    if _CHECKPOINTER is not None:
        _CHECKPOINTER.save(ordinal, name, {
            'args': args,
            'proteins': proteins,
            'all_proteins': all_proteins if all_proteins is not None else proteins,
        })
    logging.info('Stage %s %s: completed %s proteins in %.1f s',
                 ordinal, name, len(proteins), seconds)
    if _RUN_RECORD is not None:
        _RUN_RECORD.stage(ordinal, name, seconds, len(proteins))
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


def run_pipeline():
    """Runs NERVE"""
    global _CHECKPOINTER, _RUN_RECORD, _RUN_STARTED
    nerve_start = time.monotonic()
    _RUN_STARTED = nerve_start
    _STAGE_TIMINGS.clear()
    args=get_args()
    args = args._replace(working_dir=os.path.join(os.path.abspath(args.working_dir), ''))
    check_output_directory(args.working_dir, args.resume == 'True')
    os.makedirs(args.working_dir, exist_ok=True)
    configure_logging(args.working_dir, args.debug == 'True', args.resume == 'True')
    logging.info('Start NERVE 2.0%s', ' (resume attempt)' if args.resume == 'True' else '')
    logging.debug('Effective arguments: %s', args._asdict())
    source_identity = {
        'commit': os.environ.get('NERVE_SOURCE_COMMIT', 'unknown'),
        'tree': os.environ.get('NERVE_SOURCE_TREE', 'unknown'),
        'dirty': os.environ.get('NERVE_SOURCE_DIRTY', 'unknown'),
    }
    if args.resume != 'True':
        _RUN_RECORD = RunRecord(args.working_dir, args, [], source_identity, None)
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
                  'checkpoint', 'resume', 'debug'):
        signature_fields.pop(field, None)
    inputs = input_records([args.proteome1, args.proteome2])
    signature = compute_signature(
        [args.proteome1, args.proteome2], signature_fields, source_identity,
        input_info=inputs)
    _CHECKPOINTER = Checkpointer(
        args.working_dir, signature,
        enabled=args.checkpoint == 'True', resume=args.resume == 'True')
    restored = _CHECKPOINTER.restored_state()
    completed_stage = _CHECKPOINTER.bundle['ordinal'] if restored else 0
    if args.resume == 'True':
        _RUN_RECORD = RunRecord(args.working_dir, args, inputs, source_identity,
                                signature, resume=True)
    else:
        _RUN_RECORD.data.update(arguments=args._asdict(), inputs=inputs,
                                run_signature=signature)
        _RUN_RECORD.save()
    if restored:
        _RUN_RECORD.restored_through(completed_stage)
        args = restored['args']._replace(debug=args.debug)
        list_of_proteins = restored['all_proteins']
        final_proteins = restored['proteins']
        logging.info('Resuming after stage %s; earlier stages restored', completed_stage)
    logging.info('Source commit: %s', source_identity['commit'])
    logging.debug('NERVE source commit=%s tree=%s dirty=%s run_signature=%s',
                 source_identity['commit'], source_identity['tree'],
                 source_identity['dirty'], signature)
    
    # run quality control
    if completed_stage < 1:
        start=time.monotonic()
        start_stage(1, 'quality_control', 'input', completed_stage)
        original_primary, original_secondary = args.proteome1, args.proteome2
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
        _RUN_RECORD.data['qc'] = {
            'primary': qc_counts(original_primary, args.working_dir),
            'secondary': qc_counts(original_secondary, args.working_dir)
            if original_secondary else None}
        record_stage(args, 1, 'quality_control', list_of_proteins, start)
    
    # subcellular localization prediction
    start = time.monotonic()
    start_stage(2, 'localization', len(list_of_proteins), completed_stage)
    logging.debug("Subcelloc start with psortb...")
    if completed_stage < 2:
        list_of_proteins = psortb(list_of_proteins, args.working_dir, args.gram, args.proteome1)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end - start))
    record_stage(args, 2, 'localization', list_of_proteins, start)

    # Adhesin
    logging.debug("Adhesin start...")
    start=time.monotonic()
    start_stage(3, 'adhesin', len(list_of_proteins), completed_stage)
    # extract features
    if completed_stage < 3:
        list_of_proteins = extract_features(list_of_proteins, args.NERVE_dir, args.iFeature_dir, args.working_dir, args.proteome1)
        list_of_proteins = adhesin_predict(list_of_proteins, args.NERVE_dir)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    record_stage(args, 3, 'adhesin', list_of_proteins, start)
    
    # Tmhelices
    logging.debug("Tmhelices start...")
    start=time.monotonic()
    start_stage(4, 'topology', len(list_of_proteins), completed_stage)
    if completed_stage < 4:
        list_of_proteins = tmhelices(list_of_proteins, args.working_dir,
                                     threads=args.threads)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    record_stage(args, 4, 'topology', list_of_proteins, start)
    
    # Razor
    if args.razor == "True":
        logging.debug("Loop-razor start...")
        start=time.monotonic()
        start_stage(5, 'razor', len(list_of_proteins), completed_stage)
        if completed_stage < 5:
            list_of_proteins=razor(list_of_proteins, args.working_dir, args.transmemb_doms_limit, args.razlen)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end-start))
        record_stage(args, 5, 'razor', list_of_proteins, start)
    
    # Autoimmunity
    logging.debug("Autoimmunity start...")
    start=time.monotonic()
    start_stage(6, 'autoimmunity', len(list_of_proteins), completed_stage)
    if completed_stage < 6:
        list_of_proteins=autoimmunity(list_of_proteins, args.proteome1, args.working_dir, args.NERVE_dir, args.e_value,
                                      args.minlength, args.mismatch, args.substitution, threads=args.threads)
    end=time.monotonic()
    logging.debug("Done run in: {:.4f} seconds".format(end-start))
    record_stage(args, 6, 'autoimmunity', list_of_proteins, start)
    
    # Mouse immunity
    if args.mouse == "True":
        start = time.monotonic()
        start_stage(7, 'mouse', len(list_of_proteins), completed_stage)
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
        start_stage(8, 'conservation', len(list_of_proteins), completed_stage)
        logging.debug("Conservation start...")
        if completed_stage < 8:
            list_of_proteins = conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value,
                                            args.proteome1, args.proteome2,
                                            args.minlength, args.substitution, args.mismatch, threads=args.threads)
        end = time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 8, 'conservation', list_of_proteins, start)
        
    # Virulence
    #logging.debug(f'list of proteins before virulence:\n{list_of_proteins}')
    if args.virulent == "True":
        start=time.monotonic()
        start_stage(9, 'virulence', len(list_of_proteins), completed_stage)
        logging.debug("Virulence start...")
        if completed_stage < 9:
            list_of_proteins = virulent_factor_predict(list_of_proteins, args.NERVE_dir)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 9, 'virulence', list_of_proteins, start)
        
    # annotation
    if args.annotation == "True":
        start = time.monotonic()
        start_stage(10, 'annotation', len(list_of_proteins), completed_stage)
        logging.debug("Annotation start...")
        if completed_stage < 10:
            list_of_proteins = annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        end=time.monotonic()
        logging.debug("Done run in: {:.4f} seconds".format(end - start))
        record_stage(args, 10, 'annotation', list_of_proteins, start)
    
    # score
    for protein in list_of_proteins:
        protein.score = scorer(protein, args.mouse_peptides_sum_limit, args.mouse)
    
    # select
    if completed_stage < 11:
        final_proteins=list_of_proteins
    if args.select == "True":
        logging.debug("Select start...")
        start=time.monotonic()
        start_stage(11, 'select', len(list_of_proteins), completed_stage)
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
        start_stage(12, 'epitope', len(final_proteins), completed_stage)
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
    final_proteins_names = {p.id for p in final_proteins}
    discarded_proteins = [p for p in list_of_proteins if p.id not in final_proteins_names]
    output(discarded_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'))
    reason_counts = write_selection_reasons(
        list_of_proteins, args, os.path.join(args.working_dir, 'selection_reasons.csv'))
    
    nerve_end = time.monotonic()
    report = timing_report()
    if report:
        logging.info('%s attempt stage timing:\n%s',
                     'Resumed' if completed_stage else 'Fresh', report)
    logging.debug("Done: NERVE has finished its analysis in: {:.4f} seconds".format(nerve_end-nerve_start))
    logging.info('End NERVE computation successfully: %s candidates, %s discarded; '
                 'outputs in %s', len(final_proteins), len(discarded_proteins), args.working_dir)
    _CHECKPOINTER.finish()
    _RUN_RECORD.finish('success', nerve_end - nerve_start,
                       candidate_count=len(final_proteins),
                       discarded_count=len(discarded_proteins),
                       selection_reason_counts=reason_counts,
                       outputs={name: {'path': os.path.join(args.working_dir, name),
                                       'sha256': file_hash(os.path.join(args.working_dir, name))}
                                for name in ('vaccine_candidates.csv', 'discarded_proteins.csv',
                                             'selection_reasons.csv')})


def main():
    global _RUN_RECORD
    def interrupt(_signum, _frame):
        raise KeyboardInterrupt
    previous_sigterm = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, interrupt)
    try:
        run_pipeline()
    except (KeyboardInterrupt, Exception) as error:
        interrupted = isinstance(error, KeyboardInterrupt)
        logging.error('NERVE %s: %s', 'interrupted' if interrupted else 'failed', error)
        if any(handler.level == logging.DEBUG for handler in logging.getLogger().handlers):
            logging.debug('Failure details', exc_info=True)
        if _RUN_RECORD is not None:
            _RUN_RECORD.finish('interrupted' if interrupted else 'failed',
                               time.monotonic() - _RUN_STARTED, error=str(error))
        raise SystemExit(130 if interrupted else 1)
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm)

if __name__ == "__main__":
    main()
