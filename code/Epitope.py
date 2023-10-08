
"""Run epitope prediction of proteins"""

import logging
from epitopepredict import base, plotting
import math
from code.Utils import *
from code.Select import *


def epitope(final_proteins, working_dir,
            mhci_length, mhcii_length, mhci_overlap, mhcii_overlap, epitope_percentile, ep_plots, transmemb_doms_limit) -> list:
    """Module to run epitopes prediction"""
    
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    # create predictor object for mhcii
    mhcii_predictor = base.get_predictor('tepitope')
    # create predictor object for mhci
    mhci_predictor = base.get_predictor('basicmhc1')

    # set alleles to consider. Based on literature (Alessandro Sette work) we choose to use supertypes alleles
    m2alleles = base.get_preset_alleles('mhc2_supertypes')
    m1alleles = base.get_preset_alleles('mhc1_supertypes')

    # calculate score for every protein in the list
    protein_scores = []
    
    for protein in final_proteins:
        score = protein.score
        protein_scores.append(score)

        # initialize protein sequence
        #if protein.transmembrane_doms > 0:
        #    new_loop_out = protein.provide_raw_loops(transmemb_doms_limit)
        #    protein.sequence_out = new_loop_out
    

    if len(protein_scores) != 0:

        # sort scores in ascending order and find selected percentile
        sorted_scores = sorted(protein_scores)
        n = len(sorted_scores)
        percentile_index = math.ceil(n * epitope_percentile) - 1
        percentile = sorted_scores[percentile_index]
        
        for p, score in zip(final_proteins, protein_scores):
            if score >= percentile:
                # create a dir for every protein
                new_dir_path = os.path.join(working_dir, 'epitope', p.accession)
                os.makedirs(new_dir_path, exist_ok=True)
                # run predictions for MHC I and II epitopes
                # use 'threads=0' to use all available cores
                
                sequence = p.sequence
                #sequence = p.sequence_out if p.sequence_out != None else p.sequence
                mhci_epitopes = mhci_predictor.predict_sequences(sequence, alleles=m1alleles, length=mhci_length,
                                                                 verbose=False, overlap=mhci_overlap)
                                                                	    
                mhcii_epitopes = mhcii_predictor.predict_sequences(sequence, alleles=m2alleles, length=mhcii_length,
                                                                   verbose=False, overlap=mhcii_overlap)
                                                                   	
                mhci_epitopes['protein id'] = p.id
                mhcii_epitopes['protein id'] = p.id   

                # filter based on default cutoff (0.95) and default by allele method
                filtered_binders1 = mhci_epitopes.copy()
                filtered_binders2 = mhcii_epitopes.copy()
                
                if not filtered_binders1.empty:
                    filtered_binders1 = mhci_predictor.get_binders()
    
                if not filtered_binders2.empty:
                    filtered_binders2 = mhcii_predictor.get_binders()
                
                # find promiscuous binders (use names = ... to work only on filtered ones)
                pb1 = pd.DataFrame(columns=['peptide', 'pos',  'alleles', 'core', 'score'])
                if not filtered_binders1.empty:
                    pb1 = mhci_predictor.promiscuous_binders(names = filtered_binders1)
                    pb1 = pb1[pb1.alleles > 1]
                if not filtered_binders2.empty:
                    pb2 = mhcii_predictor.promiscuous_binders(names = filtered_binders2)
                    pb2 = pb2[pb2.alleles > 1]

                if not filtered_binders1.empty:
                    tmp_dic = {}
                    for allele in filtered_binders1.allele.unique():
                        seq_cov = [0 for _ in sequence]
                        for start in filtered_binders1[filtered_binders1.allele == allele].pos.to_list():
                            for pos in range(start, start+mhci_length):
                                seq_cov[pos - 1] += 1
                        tmp_dic[allele] = " ".join([str(_) for _ in seq_cov])

                    p.MHC1_binders = tmp_dic
            
                if not filtered_binders2.empty:
                    tmp_dic = {}
                    for allele in filtered_binders2.allele.unique():
                        seq_cov = [0 for _ in sequence]
                        for start in filtered_binders2[filtered_binders2.allele == allele].pos.to_list():
                            for pos in range(start, start+mhcii_length):
                                seq_cov[pos - 1] += 1
                        tmp_dic[allele] = " ".join([str(_) for _ in seq_cov])

                    p.MHC2_binders = tmp_dic

                if not pb1.empty:
                    p.MHC1_pb_binders = pb1[['peptide', 'pos',  'alleles', 'core', 'score']].reset_index(drop=True).to_csv(sep='\t')

                if not pb2.empty:
                    p.MHC2_pb_binders = pb2[['peptide', 'pos',  'alleles', 'core', 'score']].reset_index(drop=True).to_csv(sep='\t')

                # save files
                mhci_epitopes.to_csv(os.path.join(new_dir_path, 'mhci_epitopes_{}.csv'.format(p.accession)), index=False)
                mhcii_epitopes.to_csv(os.path.join(new_dir_path, 'mhcii_epitopes_{}.csv'.format(p.accession)), index=False)
                filtered_binders1.to_csv(os.path.join(new_dir_path, 'MHC1_epitopes_FILTERED_{}.csv'.format(p.accession)), index=False)
                filtered_binders2.to_csv(os.path.join(new_dir_path, 'MHC2_epitopes_FILTERED_{}.csv'.format(p.accession)), index=False)
                pb1.to_csv(os.path.join(new_dir_path, 'Promiscuous_binders_MHC1_{}.csv'.format(p.accession)), index=False)
                pb2.to_csv(os.path.join(new_dir_path, 'Promiscuous_binders_MHC2_{}.csv'.format(p.accession)), index=False)
                
                # plot binders in a sequence
                if ep_plots==True:
                
                    names_i = mhci_predictor.get_names()
                    for name in names_i:
                  
                       ax = plotting.plot_tracks([mhci_predictor], name=name)
                       if ax != None:
                          ax.figure.savefig(fname=os.path.join(new_dir_path, 'tracks_plot_pbs_MHC1_{}.png'.format(p.accession)))
                  
                   
                    names_ii = mhcii_predictor.get_names()
                    for name in names_ii:
                   
                       ax = plotting.plot_tracks([mhcii_predictor], name=name)
                       if ax != None:
                          ax.figure.savefig(fname=os.path.join(new_dir_path, 'tracks_plot_pbs_MHC2_{}.png'.format(p.accession)))
                   
                
                # plot heatmap colored by ranks
                    for name in names_i:
                   
                       ax = plotting.plot_binder_map(mhci_predictor, name=name)
                       if ax != None:
                          ax.figure.savefig(fname=os.path.join(new_dir_path, 'heatmap_pbs_MHC1_{}.png'.format(p.accession)))
                  
                
                    for name in names_ii:
                   
                       ax = plotting.plot_binder_map(mhcii_predictor, name=name)
                       if ax != None:
                          ax.figure.savefig(fname=os.path.join(new_dir_path, 'heatmap_pbs_MHC2_{}.png'.format(p.accession)))
                         
                
    return final_proteins
