#!/usr/bin/python3
"""Select module"""

import pandas as pd
from Protein import Protein
    
def select(list_of_proteins, transmemb_doms_limit,
           padlimit, mouse, mouse_peptides_sum_limit, virlimit, virulent, annotation)->list:
    """Selects suitable candidate proteins for vaccine production"""
    
    # annotations to exclude proteins
    exclusion_annotations = ['structural constituent of ribosome', 'DNA binding', 'DNA-binding transcription factor activity',
                             'transcription regulator activity', 'rRNA binding', 'RNA binding',
                             'aminoacyl-tRNA ligase activity', 'sequence-specific DNA binding', 
                             'catalytic activity, acting on a tRNA', 'catalytic activity, acting on RNA',
                             'tyrosine-tRNA ligase activity', 'aminoacyl-tRNA editing activity',
                             'translation factor activity, RNA binding', 'translation regulator activity',
                             'translation regulator activity, nucleic acid binding', 
                             'translation elongation factor activity', 'catalytic activity, acting on DNA'
                            ]
    
    final_list = []
    for protein in list_of_proteins:
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
        
        annotation_flag = "False"
        if annotation == "True":
            for annot in exclusion_annotations:
                if annot in str(protein.annotations): 
                    annotation_flag = "True"
        if annotation_flag == "True": continue

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
                 (round(scorer(protein, mouse_peptides_sum_limit, mouse), 4)),
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
                 str("".join([str(round(protein.mouse_peptides_sum,4)) if protein.mouse_peptides_sum!=None else "0"])),
                 str("".join([str(protein.annotations) if protein.annotations!=None else ""])),
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)))), 
                 str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)))),  
                 str(protein.sequence),
                 str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor!=None else ""])),
                 str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""])), # should be shown anyways
                 str("".join([str(protein.HLA_A_01_01) if protein.HLA_A_01_01!=None else ""])),
                 str("".join([str(protein.HLA_A_02_01) if protein.HLA_A_02_01!=None else ""])),
                 str("".join([str(protein.HLA_A_03_01) if protein.HLA_A_03_01!=None else ""])),
                 str("".join([str(protein.HLA_A_24_02) if protein.HLA_A_24_02!=None else ""])),
                 str("".join([str(protein.HLA_B_07_02) if protein.HLA_B_07_02!=None else ""])),
                 str("".join([str(protein.HLA_B_44_03) if protein.HLA_B_44_03!=None else ""])),
                 str("".join([str(protein.HLA_A_01_01) if protein.HLA_A_01_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_01_01) if protein.HLA_DRB1_01_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_03_01) if protein.HLA_DRB1_03_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_04_01) if protein.HLA_DRB1_04_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_07_01) if protein.HLA_DRB1_08_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_11_01) if protein.HLA_DRB1_11_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_13_01) if protein.HLA_DRB1_13_01!=None else ""])),
                 str("".join([str(protein.HLA_DRB1_15_01) if protein.HLA_DRB1_15_01!=None else ""])),
                 str("".join([str(protein.pb1) if protein.pb1!=None else ""])),
                 str("".join([str(protein.pb2) if protein.pb2!=None else ""]))
                 
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
                    'MHC1_HLA_A_01_01',
                    'MHC1_HLA_A_02_01',
                    'MHC1_HLA_A_03_01',
                    'MHC1_HLA_A_24_02',
                    'MHC1_HLA_B_07_02',
                    'MHC1_HLA_B_44_03',
                    'MHC2_HLA_DRB1_01_01',
                    'MHC2_HLA_DRB1_03_01',
                    'MHC2_HLA_DRB1_04_01',
                    'MHC2_HLA_DRB1_07_01',
                    'MHC2_HLA_DRB1_08_01',
                    'MHC2_HLA_DRB1_11_01',
                    'MHC2_HLA_DRB1_13_01',
                    'MHC2_HLA_DRB1_15_01',
                    'promiscuous_binders_MHC1',
                    'promiscuous_binders_MHC2'
                     ]
                )
    df = df.sort_values(by = 'score', ascending = False)
    df.to_csv(outfile, index = False)
