#!/usr/bin/python3

import requests, os
from typing import NamedTuple
import argparse

class Args(NamedTuple):
    '''Command-line arguments'''
    proteome_id:str
    output_dir:str

def get_args() -> Args:
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser(
        description='Downloads proteome and saves in fasta format',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-proteome_id',
                        metavar='--proteome-id', 
                        help='Uniprot proteome accession code',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-output_dir',
                        metavar='--output-dir', 
                        help='Output directory',
                        type=str,
                        required=False,
                        )
    args = parser.parse_args()
    return Args(args.proteome_id, args.output_dir)

def main():
    args = get_args()
    proteome_downloader(args.proteome_id, args.output_dir)

def proteome_downloader(proteome_id, output_dir = None, format = "fasta") -> None:
  """Downloads proteome from uniprot database
  param: proteome_id: uniprot unique proteome id"""
  url = f'https://www.uniprot.org/uniprot/?query=proteome:{proteome_id}&format={format}' 
  filename = f'{proteome_id}_proteome.fasta'
  output_dir = output_dir if output_dir != None else os.getcwd()
  try:
     response = requests.get(url, stream = True)
     text_file = open(os.path.join(output_dir, filename), 'wb')
     for chunk in response.iter_content(chunk_size=1024):
        text_file.write(chunk)
     text_file.close()
     print(f'Proteome {proteome_id} downloaded succesfully')
     #output = open(os.path.join(output_dir, filename), 'r').readlines()
     #return output
  except Exception as e:
    print(f'An error has occured:\n{e}')
    return None
    
if __name__ == "__main__":
    main()
