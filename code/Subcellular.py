#!/usr/bin/python3
"""Runs subcellular localization predictions with psortb"""

import os, logging, urllib, json
from operator import attrgetter

class Localization:
        """class to store and handle protein subcellular localizations"""
        def __init__(self, localization, reliability):
            self.localization = str(localization)
            self.reliability = float(reliability)

def psortb(list_of_proteins, working_dir, gram, proteome1) -> list:
    """Psortb docker-based api request sender and  result parser"""
    # define logging file
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)
    
    url = "http://psortb:8080/execute"
    # read proteome1 file and pass it to json request
    infile = open(proteome1, 'r')
    proteome1 = "".join(infile)
    body = {'gram':gram,'seq':proteome1}
    with open(os.path.join(working_dir, 'payload.json'), 'w') as f:
        json.dump(body, f)
    req = urllib.request.Request(url)
    req.add_header('Content-Type', 'application/json; charset=utf-8')
    jsondata = json.dumps(body)
    jsondataasbytes = jsondata.encode('utf-8')   # needs to be bytes
    req.add_header('Content-Length', len(jsondataasbytes))
    logging.debug("Sending request to psortb")
    response = urllib.request.urlopen(req, jsondataasbytes)
    data_json = json.loads(response.read())
    logging.debug(f'Psortb stdout:\n{data_json["stdout"]}\nPsortb stderr:\n{data_json["stderr"]}')
    logging.debug('Parsing psortb output')
    for entry in data_json['result'].split('-------------------------------------------------------------------------------\n\n'):
        split = entry.split('\n')
        id_ = split[0][split[0].find('SeqID: ')+len('SeqID: '):].strip() # protein id
        # extract sublocalization predictions:
        if 'Localization Scores:' in entry:
            # get part of the output that define the predictions
            predictions = entry[entry.find('Localization Scores:') + len('Localization Scores:'):entry.find('Final')]
            # get list containing each prediction and its score
            predictions = (predictions.strip().split('\n'))
            localizations = [Localization(element.split()[0], element.split()[1]) for element in predictions]
            # sort localizations
            localizations = sorted(localizations, key = attrgetter('reliability'), reverse = True)
            # set unknown prediction
            if localizations[0].reliability <= 3.:
                localizations = [Localization('Unknown', 0)]
            for p in list_of_proteins:
                if id_ in p.id:
                    p.localization = localizations
    # save psortb raw predictions
    #df.to_csv(os.path.join(working_dir, 'psortb_predictions.csv'))
    return list_of_proteins
