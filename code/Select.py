#!/usr/bin/python3
"""Select module"""

import csv
import os
from collections import Counter

import pandas as pd

from code.Protein import Protein


def selection_reason(protein, transmemb_doms_limit, padlimit, mouse,
                     mouse_peptides_sum_limit, virlimit, virulent, razor):
    """First failed selection rule, or None."""
    location = protein.localization[0]
    if location.localization == 'Cytoplasmic':
        return 'cytoplasmic', location.localization, 'non-cytoplasmic'
    if location.reliability < 7.49:
        return 'localization_score', location.reliability, 7.49
    if virulent == 'False' and protein.p_ad < padlimit:
        return 'adhesin_probability', protein.p_ad, padlimit
    if virulent == 'True' and protein.p_ad < padlimit and protein.p_vir < virlimit:
        return ('adhesin_and_virulence_probability',
                f'{protein.p_ad},{protein.p_vir}', f'{padlimit},{virlimit}')
    if protein.transmembrane_doms >= transmemb_doms_limit:
        if razor == 'False' or protein.original_sequence_if_razor is None:
            return 'transmembrane_domains', protein.transmembrane_doms, transmemb_doms_limit
    if protein.sapiens_peptides_sum > .15:
        return 'human_peptide_sum', protein.sapiens_peptides_sum, .15
    if protein.list_of_peptides_from_comparison_with_mhcpep_sapiens:
        return 'human_mhcpep_match', len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens), 0
    if mouse == 'True':
        if protein.mouse_peptides_sum > mouse_peptides_sum_limit:
            return 'mouse_peptide_sum', protein.mouse_peptides_sum, mouse_peptides_sum_limit
        if protein.list_of_peptides_from_comparison_with_mhcpep_mouse:
            return 'mouse_mhcpep_match', len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse), 0
    return None
    
def select(list_of_proteins, transmemb_doms_limit,
           padlimit, mouse, mouse_peptides_sum_limit, virlimit, virulent, razor)->list:
    """Selects suitable candidate proteins for vaccine production"""
        
    final_list = []
    for protein in list_of_proteins:
        if selection_reason(protein, transmemb_doms_limit, padlimit, mouse,
                            mouse_peptides_sum_limit, virlimit, virulent, razor) is None:
            final_list.append(protein)
    return final_list


def write_selection_reasons(proteins, args, path):
    """Write selection decisions."""
    counts = Counter()
    path = os.fspath(path)
    temporary = path + '.tmp'
    try:
        with open(temporary, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['accession', 'selected', 'first_rule', 'observed', 'cutoff'])
            for protein in proteins:
                reason = (selection_reason(protein, args.transmemb_doms_limit,
                                           args.padlimit, args.mouse,
                                           args.mouse_peptides_sum_limit,
                                           args.virlimit, args.virulent, args.razor)
                          if args.select == 'True' else None)
                rule = reason[0] if reason else ('selected' if args.select == 'True' else 'selection_disabled')
                counts[rule] += 1
                writer.writerow([protein.accession, reason is None,
                                 rule, reason[1] if reason else '',
                                 reason[2] if reason else ''])
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)
    return dict(counts)


def scorer(protein:Protein, mouse_peptides_sum_limit:float, mouse:str) -> float:
    """Provides a score for protein candidates"""
    
    if mouse == "True":
        score = (protein.p_ad + (protein.p_vir if protein.p_vir != None else 1 ) +\
             ((protein.localization[0].reliability / 10) if protein.localization[0].reliability > 7.49 else 0) +\
             (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) +\
             (1 - (protein.sapiens_peptides_sum / .15)) + (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)) +\
             (1 - (protein.mouse_peptides_sum / mouse_peptides_sum_limit))) / 7
    if mouse != "True":
        score = (protein.p_ad + (protein.p_vir if protein.p_vir != None else 1 ) +\
             ((protein.localization[0].reliability / 10) if protein.localization[0].reliability > 7.49 else 0) +\
             (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) +\
             (1 - (protein.sapiens_peptides_sum / .15))) / 5
    return score

def output(list_of_proteins:list, outfile):
    """Produces output .csv table"""
    df = pd.DataFrame([[str(getattr(protein, 'original_id', protein.id)),
                 str("".join([str(protein.accession) if protein.accession!=None else ""])),
                 (round(protein.score, 4)),
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
                 # Remove duplicates while preserving order.
                 str(", ".join(dict.fromkeys(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens))),
                 str(", ".join(dict.fromkeys(protein.list_of_peptides_from_comparison_with_mhcpep_mouse))),
                 str(protein.sequence),
                 str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor!=None else ""])),
                 #str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""])),
                 str(str(protein.tmhmm_seq)),
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
                    'transmembrane_domains',
                    'localization',
                    'localization score',
                    'virulence_probability',
                    'adhesin_probability',
                    'conservation_score',
                    'shared_human_peptides',
                    'shared_mouse_peptides',
                    'shared_conserved_proteome_peptides',
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
                    'MHC1_promiscuous_binders',
                    'MHC2_promiscuous_binders',
                   
                     ]
                )
    df = df.sort_values(by = 'score', ascending = False)
    df.to_csv(outfile, index = False)
