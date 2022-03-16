# #!/usr/bin/python3

import os
from urllib.parse import urlencode
import pycurl
from io import BytesIO 
from bs4 import BeautifulSoup as bs
from typing import NamedTuple
from operator import attrgetter
import pandas as pd
import argparse

class Args(NamedTuple):
    '''Command-line arguments'''
    path_to_fastas:str
    gram_type:str

class DeepFriEntry(NamedTuple):
    '''Handles DeepFri entry'''
    id_:str
    score:float
    GO:str

def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description='Downloads cello subcellular localization predictions. For each protein the prediction with the highest score (and the score itself) are reported.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-path_to_fastas',
                        metavar='--path-to-fastas', 
                        help='Path where sequences to evaluate are stored',
                        type=open,
                        required=True,
                        )
    parser.add_argument('-gram_type',
                        metavar='--gram-type', 
                        help='Negative/Positive',
                        type=str,
                        required=True,
                        )
    args = parser.parse_args()
    return Args(args.path_to_fastas, args.gram_type)

def main():
  args = get_args() 
  cello_text_output = cello_scraper(args.gram_type, args.path_to_fastas)
  df = cello_output_parser(cello_text_output)
  df.to_excel('./cello_output.xlsx', index=False)
  print(df)

def cello_scraper(gramtype:str, infile:str) -> str:
  """Downloads protein localization predictions from CELLO (http://cello.life.nctu.edu.tw/) webpage
  param: gramtype: 'Negative' or 'Positive'
  param: inifle: path to input fasta file"""
  # set connection time and operation time to 1 hour
  curlConnector = pycurl.Curl()
  curlConnector.setopt(pycurl.CONNECTTIMEOUT, 3600)
  curlConnector.setopt(pycurl.TIMEOUT, 3600)
  # set variables
  gramtype = 'pro' if gramtype == 'Positive' else 'gramp' if gramtype == 'Negative' else None
  infile = infile.readlines() # infile is already "open" at this point
  proteome = "".join(infile)
  # pycurl recipe from http://pycurl.io/docs/latest/quickstart.html:
  b_obj = BytesIO() 
  crl = pycurl.Curl()
  crl.setopt(crl.URL, 'http://cello.life.nctu.edu.tw/cgi/main.cgi') # <---------------------------------
  # Write bytes that are utf-8 encoded
  crl.setopt(crl.WRITEDATA, b_obj)
  data = {'species':gramtype, 'seqtype':'prot', 'fasta':proteome, 'Submit':'Submit'} # <----------------
  pf = urlencode(data)
  # Sets request method to POST,
  # Content-Type header to application/x-www-form-urlencoded
  # and data to send in request body.
  crl.setopt(crl.POSTFIELDS, pf)
  crl.perform()
  crl.close()
  get_body = b_obj.getvalue().decode('utf8')
  # convert HTML to nice text
  soup = bs(get_body)  
  soup = soup.get_text()
  return soup

def cello_output_parser(cello_text_output:str)->pd.DataFrame():
  """Parses cello output"""
  # define class to handle predictions
  Localization = namedtuple('Localization',['localization','reliability'])
  outlist=[]
  for line in (cello_text_output.split('*********************************************************************************')):
    if 'Documentation' not in line:
       localizations = []
       seq_id = line[line.find('SeqID')+7:line.find('Analysis')]
       prediction = line[line.find('Prediction')+15:].split('\xa0\xa0\xa0\xa0')
       # collect only significative predictions (*-marked)
       for i in range(0, len(prediction),2):
         if '*' in prediction[i+1]:
            localizations.append(Localization(prediction[i][1:].split('\t')[0], float(prediction[i+1][:4])))
       # sort predictions based on reliability score
       sorted_predictions = sorted(localizations, key=attrgetter('reliability'), reverse=True)
       outlist.append([seq_id, sorted_predictions])
  df = pd.DataFrame(outlist, index=range(len(outlist)))
  return df
    
if __name__ == "__main__":
    main() 
