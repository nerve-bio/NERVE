"""Parallel epitope prediction with one reusable predictor pair per worker."""
import math, multiprocessing, os
from concurrent.futures import ProcessPoolExecutor
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ.setdefault(_name, '1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from epitopepredict import base, plotting
from code.Utils import epitope_output_identifier

_PREDICTORS = {}

def available_cpu_count():
    """CPU count constrained by affinity and Docker cgroup quotas."""
    counts = [os.cpu_count() or 1]
    if hasattr(os, 'sched_getaffinity'):
        counts.append(len(os.sched_getaffinity(0)))
    try:
        if os.path.isfile('/sys/fs/cgroup/cpu.max'):
            quota, period = open('/sys/fs/cgroup/cpu.max').read().split()[:2]
            if quota != 'max': counts.append(max(1, int(float(quota) / float(period))))
        elif os.path.isfile('/sys/fs/cgroup/cpu/cpu.cfs_quota_us'):
            quota = int(open('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read())
            period = int(open('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read())
            if quota > 0: counts.append(max(1, int(quota / period)))
    except (OSError, ValueError, ZeroDivisionError):
        pass
    return max(1, min(counts))

def resolve_epitope_workers(requested, task_count):
    if requested is None or requested == 0: requested = available_cpu_count()
    if requested < 0: raise ValueError('epitope_workers must be zero (auto) or positive')
    return max(1, min(requested, max(1, task_count)))

def _init_worker():
    _PREDICTORS['mhcii'] = base.get_predictor('tepitope')
    _PREDICTORS['mhci'] = base.get_predictor('basicmhc1')
    _PREDICTORS['m2'] = base.get_preset_alleles('mhc2_supertypes')
    _PREDICTORS['m1'] = base.get_preset_alleles('mhc1_supertypes')

def _coverage(binders, sequence, length):
    result = {}
    for allele in binders.allele.unique():
        coverage = [0 for _ in sequence]
        for start in binders[binders.allele == allele].pos.to_list():
            for position in range(start, min(start + length, len(sequence) + 1)):
                coverage[position - 1] += 1
        result[allele] = ' '.join(map(str, coverage))
    return result

def _save_plots(predictor, output_dir, prot_name, kind):
    for name in predictor.get_names():
        ax = plotting.plot_tracks([predictor], name=name)
        if ax is not None:
            ax.figure.savefig(os.path.join(output_dir, f'tracks_plot_pbs_{kind}_{prot_name}.png'))
            plt.close(ax.figure)
        ax = plotting.plot_binder_map(predictor, name=name)
        if ax is not None:
            ax.figure.savefig(os.path.join(output_dir, f'heatmap_pbs_{kind}_{prot_name}.png'))
            plt.close(ax.figure)

def _predict_candidate(task):
    p, working_dir, m1len, m2len, m1overlap, m2overlap, ep_plots = task
    if not _PREDICTORS: _init_worker()
    m1, m2 = _PREDICTORS['mhci'], _PREDICTORS['mhcii']
    sequence = p.sequence
    name = epitope_output_identifier(p)
    output_dir = os.path.join(working_dir, 'epitope', name)
    os.makedirs(output_dir, exist_ok=True)
    mhci = m1.predict_sequences(sequence, alleles=_PREDICTORS['m1'], length=m1len, verbose=False, overlap=m1overlap)
    mhcii = m2.predict_sequences(sequence, alleles=_PREDICTORS['m2'], length=m2len, verbose=False, overlap=m2overlap)
    mhci['protein id'], mhcii['protein id'] = p.id, p.id
    binders1 = m1.get_binders() if not mhci.empty else mhci.copy()
    binders2 = m2.get_binders() if not mhcii.empty else mhcii.copy()
    columns = ['peptide','pos','alleles','core','score']
    pb1, pb2 = pd.DataFrame(columns=columns), pd.DataFrame(columns=columns)
    if not binders1.empty:
        pb1 = m1.promiscuous_binders(names=binders1); pb1 = pb1[pb1.alleles > 1]
        p.MHC1_binders = _coverage(binders1, sequence, m1len)
    if not binders2.empty:
        pb2 = m2.promiscuous_binders(names=binders2); pb2 = pb2[pb2.alleles > 1]
        p.MHC2_binders = _coverage(binders2, sequence, m2len)
    if not pb1.empty: p.MHC1_pb_binders = pb1[columns].reset_index(drop=True).to_csv(sep='\t')
    if not pb2.empty: p.MHC2_pb_binders = pb2[columns].reset_index(drop=True).to_csv(sep='\t')
    for frame, filename in ((mhci,f'mhci_epitopes_{name}.csv'),(mhcii,f'mhcii_epitopes_{name}.csv'),
        (binders1,f'MHC1_epitopes_FILTERED_{name}.csv'),(binders2,f'MHC2_epitopes_FILTERED_{name}.csv'),
        (pb1,f'Promiscuous_binders_MHC1_{name}.csv'),(pb2,f'Promiscuous_binders_MHC2_{name}.csv')):
        frame.to_csv(os.path.join(output_dir, filename), index=False)
    if ep_plots == 'True':
        _save_plots(m1, output_dir, name, 'MHC1'); _save_plots(m2, output_dir, name, 'MHC2')
    return p, os.getpid()

def epitope(final_proteins, working_dir, mhci_length, mhcii_length, mhci_overlap,
            mhcii_overlap, epitope_percentile, ep_plots, transmemb_doms_limit,
            epitope_workers=0):
    del transmemb_doms_limit
    if not final_proteins: return final_proteins
    scores = [p.score for p in final_proteins]
    cutoff = sorted(scores)[math.ceil(len(scores) * epitope_percentile) - 1]
    indices = [i for i,(p,s) in enumerate(zip(final_proteins,scores))
               if s >= cutoff and len(p.sequence) >= max(mhci_length,mhcii_length)]
    tasks = [(final_proteins[i],working_dir,mhci_length,mhcii_length,mhci_overlap,mhcii_overlap,ep_plots) for i in indices]
    workers = resolve_epitope_workers(epitope_workers, len(tasks))
    if workers == 1:
        _init_worker(); results = [_predict_candidate(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn'), initializer=_init_worker) as executor:
            results = list(executor.map(_predict_candidate, tasks))
    for index,(protein,_) in zip(indices,results): final_proteins[index] = protein
    if results:
        print(f'epitope workers: requested={epitope_workers}, resolved={workers}, processes_used={len(set(pid for _,pid in results))}')
    return final_proteins
