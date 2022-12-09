#!/usr/bin/python3
"""Select module"""

import pandas as pd

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
        # proteins with Unknown localization have score==0
        #if (protein.localization[0].localization == "Unknown") and (protein.p_ad < padlimit): continue
        if mouse==True:
            if protein.mouse_peptides_sum > mouse_peptides_sum_limit: continue 
            if len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse) >= 1: continue 
        if virulent==True:
            if protein.p_vir < virlimit: continue
        final_list.append(protein)
    return final_list

def output(list_of_proteins, outfile):
    """Produces output .csv table"""
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
                    'tmhmm_seq'
                     ]
                ).to_csv(outfile)
