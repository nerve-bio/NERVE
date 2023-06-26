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
        if transmem_doms_limit = 0:
            if protein.transmembrane_doms >= transmem_doms_limit:
            
                new_loop = protein.provide_raw_loops(transmem_doms_limit)
        else
            if transmem_doms_limit > 0:
                if protein.transmembrane_doms >= transmem_doms_limit:
                    new_loop = max(p.provide_raw_loops_std(), key = lambda k: len(k))
            
            
            
            if len(new_loop) > min_loop_length:
                logging.debug(f'Substituting {str(protein.id)} sequence with its outer loops')
                protein.original_sequence_if_razor = protein.sequence
                protein.sequence = new_loop
                protein.razored = True
            else:
                logging.debug(f"No replacement found for {str(protein.id)}")
                protein.razored = False
        
    return list_of_proteins
