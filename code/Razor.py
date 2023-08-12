#!/usr/local/bin/python
"""Runs loop razor"""

import logging, os

def razor(list_of_proteins, working_dir, transmemb_doms_limit, razlen)->list:
    "Runs razor module"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
                        
    #logging.debug("Warning: razor uses X as an exclusive symbol to split the final protein. Check if X is used inside the protein sequence!")
    for protein in list_of_proteins:
        if transmemb_doms_limit:
            if protein.transmembrane_doms > transmemb_doms_limit:
            
                new_loop = protein.provide_raw_loops(transmemb_doms_limit)
            
            
            
            
            if len(new_loop) > razlen:
                logging.debug(f'Substituting {str(protein.id)} sequence with its outer loops')
                protein.original_sequence_if_razor = protein.sequence
                protein.sequence = new_loop
                protein.razored = True
            else:
                logging.debug(f"No replacement found for {str(protein.id)}")
                protein.razored = False
        
    return list_of_proteins
