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
    for p in list_of_proteins:
        if p.transmembrane_doms >= transmemb_doms_limit:
            longest_loop = max(p.provide_raw_loops(), key = lambda k: len(k))
            if len(longest_loop) > razlen:
                logging.debug(f"Substituting {str(p.id)} sequence with its longest loop.")
                p.original_sequence_if_razor = p.sequence
                p.sequence = longest_loop
    return list_of_proteins
