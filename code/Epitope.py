
"""Run epitope prediction of proteins"""


from epitopepredict import base, sequtils, analysis, plotting
import math
from utils import *
from Select import *


def epitope(final_proteins, autoimmunity, mouse, mouse_peptides_sum_limit, antigen, working_dir,
            mhci_length, mhcii_length, mhci_overlap, mhcii_overlap, epitope_percentile) -> list:
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
    #data_i = []
    #data_ii = []
    for protein in final_proteins:
        score = scorer(protein, mouse_peptides_sum_limit, mouse, autoimmunity, antigen)
        protein_scores.append(score)

    if len(protein_scores) != 0:

        # sort scores in ascending order and find selected percentile
        sorted_scores = sorted(protein_scores)
        n = len(sorted_scores)
        percentile_index = math.ceil(n * epitope_percentile) - 1
        percentile = sorted_scores[percentile_index]
	
	
        for p, score in zip(final_proteins, protein_scores):
            if score >= epitope_percentile:
                # create a dir for every protein
                new_dir_path = working_dir+'{}/'.format(p.accession)
                if os.path.isfile(new_dir_path):
                    os.removedirs(new_dir_path)
                else:
                    os.makedirs(new_dir_path)
                # run predictions for MHC I and II epitopes
                # use 'threads=0' to use all available cores
                mhci_epitopes = mhci_predictor.predict_sequences(p.sequence, alleles=m1alleles, length=mhci_length,
                                                                 verbose=False, overlap=mhci_overlap)
                                                                	
                                                                 
                mhcii_epitopes = mhcii_predictor.predict_sequences(p.sequence, alleles=m2alleles, length=mhcii_length,
                                                                   verbose=False, overlap=mhcii_overlap)
                                                                   	
                mhci_epitopes['protein id'] = p.id
                # data_i.append(mhci_epitopes)
                mhcii_epitopes['protein id'] = p.id
                # data_ii.append(mhcii_epitopes)
                
                mhci_epitopes.to_csv(new_dir_path+'mhci_epitopes_{}.csv'.format(p.accession), index=False)
                mhcii_epitopes.to_csv(new_dir_path+'mhcii_epitopes_{}.csv'.format(p.accession), index=False)
                
                p.mhci_epitopes = 'Epitopes acquired. See mhci_epitopes.csv'
                p.mhcii_epitopes = 'Epitopes acquired. See mhcii_epitopes.csv'
                
                # promiscuous binders mhc1
                results_mhc1_raw = base.results_from_csv(path=new_dir_path+'mhci_epitopes_{}.csv'.format(p.accession))
                ###
                
                score_threshold = results_mhc1_raw['score'].quantile(0.95)
                filtered_binders1 = results_mhc1_raw.loc[results_mhc1_raw['score'] >= score_threshold]
                ###
                
                #filtered_binders1 = mhci_predictor.get_binders(names=results_mhc1_raw, cutoff=0.95)  #non funziona, vedi sopra
                
                #save filtered binders
                filtered_binders1.to_csv(new_dir_path+'MHC1_epitopes_FILTERED{}.csv'.format(p.accession), index=False) #####
                # find promiscuous binders
                pb1 = mhci_predictor.promiscuous_binders(cutoff=.95, cutoff_method='score')  #cutoff=.95, cutoff_method='score'
                # save pbs
                pb1.to_csv(new_dir_path+'Promiscuous_binders_MHC1_{}.csv'.format(p.accession), index=False)
                
                # promiscuous binders mhc2
                results_mhc2_raw = base.results_from_csv(path=new_dir_path+'mhcii_epitopes_{}.csv'.format(p.accession))
                filtered_binders2 = mhcii_predictor.get_binders(names=results_mhc2_raw, cutoff=0.95)
                # save filtered binders
                filtered_binders2.to_csv(new_dir_path+'MHC2_epitopes_FILTERED{}.csv'.format(p.accession), index=False)
                # find promiscuous binders
                pb2 = mhcii_predictor.promiscuous_binders(cutoff=0.95)
                # save pbs
                pb2.to_csv(new_dir_path+'Promiscuous_binders_MHC2_{}.csv'.format(p.accession), index=False)
                
                # plot binders in a sequence
                names_i = mhci_predictor.get_names()
                for name in names_i:
                   ax = plotting.plot_tracks([mhci_predictor], name=name)
                   ax.figure.savefig(fname=new_dir_path+'tracks_plot_pbs_MHC1_{}.png'.format(p.accession))
                   
                names_ii = mhcii_predictor.get_names()
                for name in names_ii:
                   ax = plotting.plot_tracks([mhcii_predictor], name=name)
                   ax.figure.savefig(fname=new_dir_path+'tracks_plot_pbs_MHC2_{}.png'.format(p.accession))
                
                # plot heatmap colored by ranks
                for name in names_i:
                   ax = plotting.plot_binder_map(mhci_predictor, name=name)
                   ax.figure.savefig(fname=new_dir_path+'heatmap_pbs_MHC1_{}.png'.format(p.accession))
                
                for name in names_ii:
                   ax = plotting.plot_binder_map(mhcii_predictor, name=name)
                   ax.figure.savefig(fname=new_dir_path+'heatmap_pbs_MHC2_{}.png'.format(p.accession))
                                   
                
    return final_proteins
