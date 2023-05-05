#!/usr/local/bin/python
"""Runs prediction on protein topology with TMHMM"""

import logging, tmhmm, os

def tmhelices(list_of_proteins, working_dir)->list:
    "Runs Tmhelices"
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    for p in list_of_proteins:
        annotation, _ = tmhmm.predict(p.sequence)
        p.tmhmm_seq = annotation
        transmembrane_domains = 0
        for i in range(len(annotation)-1):
            if (annotation[i] == 'i' or annotation[i] == 'o') and annotation[i+1] == 'M':
                transmembrane_domains += 1
        p.transmembrane_doms = transmembrane_domains
    return list_of_proteins
