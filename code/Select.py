#!/usr/bin/python3
"""Select module"""

import pandas as pd
from code.Protein import Protein
    
def select(list_of_proteins, transmemb_doms_limit,
           padlimit, mouse, mouse_peptides_sum_limit, virlimit, virulent, razor)->list:
    """Selects suitable candidate proteins for vaccine production"""
        
    final_list = []
    for protein in list_of_proteins:

        if protein.localization[0].localization == "Cytoplasmic" or protein.localization[0].reliability < 7.49: continue
        if virulent == 'False':
            if protein.p_ad < padlimit: continue
        if virulent == 'True':
            if protein.p_ad < padlimit and protein.p_vir < virlimit: continue
        if razor == 'True':
            if (protein.transmembrane_doms >= transmemb_doms_limit) and (protein.original_sequence_if_razor is None): continue
        if razor == 'False':
            if protein.transmembrane_doms >= transmemb_doms_limit: continue

        # exclude cytoplasmatic proteins if low PAD or VIR
        if virulent == "True":
            if (protein.localization[0].localization == "Cytoplasmic" and (protein.p_ad < padlimit and protein.p_vir < virlimit)): continue 
        if virulent != "True":
            if protein.localization[0].localization == "Cytoplasmic" and protein.p_ad < padlimit: continue 
        if protein.localization[0].localization == "Cytoplasmic" and protein.localization[0].reliability >= 7.49: continue
        
        # exlude low fidelity localization prediction proteins if low PAD or VIR
        if virulent == "True":
            if protein.localization[0].reliability < 7.49 and (protein.p_vir < virlimit and protein.p_ad < padlimit): continue
        if virulent != "True":
            if protein.localization[0].reliability < 7.49 and protein.p_ad < padlimit: continue
        
        if (protein.transmembrane_doms >= transmemb_doms_limit) and (protein.original_sequence_if_razor is None): continue

        if protein.sapiens_peptides_sum > .15: continue
        if len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens) >= 1: continue
        if mouse == "True":
            if protein.mouse_peptides_sum > mouse_peptides_sum_limit: continue 
            if len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse) >= 1: continue

        final_list.append(protein)
    return final_list 

def scorer(protein:Protein, mouse_peptides_sum_limit:float, mouse:str) -> float:
    """Provides a score for protein candidates"""
    
    if mouse == "True":
        score = (protein.p_ad + (protein.p_vir if protein.p_vir != None else 0 ) +\
             ((protein.localization[0].reliability / 10) if protein.localization[0].reliability > 7.49 else 0) +\
             (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) +\
             (1 - (protein.sapiens_peptides_sum / .15)) + (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)) +\
             (1 - (protein.mouse_peptides_sum / mouse_peptides_sum_limit))) / 7
    if mouse != "True":
        score = (protein.p_ad + (protein.p_vir if protein.p_vir != None else 0 ) +\
             ((protein.localization[0].reliability / 10) if protein.localization[0].reliability > 7.49 else 0) +\
             (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) +\
             (1 - (protein.sapiens_peptides_sum / .15))) / 7
    return score

def output(list_of_proteins:list, outfile, mouse_peptides_sum_limit:float, mouse:str):
    """Produces output .csv table"""
    df = pd.DataFrame([[str(protein.id),
                 str("".join([str(protein.accession) if protein.accession!=None else ""])),
                 (round(protein.score, 4)),
                 str(protein.length),
                 str(protein.transmembrane_doms),
                 str(protein.localization[0].localization),
                 str(protein.localization[0].reliability),
                 str(", ".join([str(element) for element in protein.localization])),
                 str("".join([str(round(protein.p_vir,4)) if protein.p_vir!=None else ""])),
                 str("".join([str(round(protein.p_ad, 4)) if protein.p_ad!=None else ""])),
                 str("".join([str(round(protein.conservation_score, 4)) if protein.conservation_score!=None else ""])),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_human_peps if len(protein.list_of_shared_human_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_mouse_peps if len(protein.list_of_shared_mouse_peps)>0])))),
                 str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_conserv_proteome_peps if len(protein.list_of_shared_conserv_proteome_peps)>0])))),
                 str("".join([str(round(protein.sapiens_peptides_sum,4)) if protein.sapiens_peptides_sum!=None else "0"])),
           
                 str("".join([str(protein.annotations) if protein.annotations!=None else ""])),
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)))), 
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)))),  
                 str(protein.sequence),
                 str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor!=None else ""])),
                 str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""])), # should be shown anyways
          
                 str("".join([str(protein.MHC1_binders) if str(protein.MHC1_binders) != None else ''])),
                 str("".join([str(protein.MHC2_binders) if str(protein.MHC2_binders) != None else ''])),
                 str("".join([str(protein.MHC1_pb_binders) if str(protein.MHC1_pb_binders) != None else ''])),
                 str("".join([str(protein.MHC2_pb_binders) if str(protein.MHC2_pb_binders) != None else ''])),
               
                 
                 
                 ] for protein in list_of_proteins
                ], 
                columns= ['id',
                    'uniprot_accession_code',
                    'score',
                    'length',
                    'transmembrane_doms',
                    'localization',
                    'localization score',
                    'virulence_probability',
                    'adhesin_probability',
                    'conservation_score',
                    'shared_human_peps',
                    'shared_mouse_peps',
                    'shared_conserv_proteome_peps',
                    'human_peptides_sum',
                    'mouse_peptides_sum',
                    'annotations',
                    'list_of_peptides_from_comparison_with_mhcpep_sapiens',
                    'list_of_peptides_from_comparison_with_mhcpep_mouse',
                    'sequence',
                    'original_sequence_if_razor',
                    'tmhmm_seq',
                    'MHC1_binders',
                    'MHC2_binders',
                    'MHC1_pb_binders',
                    'MHC2_pb_binders',
                   
                     ]
                )
    df = df.sort_values(by = 'score', ascending = False)
    df.to_csv(outfile, index = False)
