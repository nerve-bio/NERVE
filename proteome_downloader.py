#!/usr/bin/python3

import requests, os
from typing import NamedTuple
import argparse

class Args(NamedTuple):
    '''Command-line arguments'''
    proteome_id:str
    output_dir:str
    out_filename:str

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
                        default=os.getcwd(),
                        required=False,
                        )
    parser.add_argument('-out_filename',
                        metavar='--output-filename', 
                        help='Output filename',
                        default='input_proteome.fasta',
                        type=str,
                        required=False,
                        )
    args = parser.parse_args()
    return Args(args.proteome_id, args.output_dir, args.out_filename)

def main():
    args = get_args()
    proteome_downloader(args.proteome_id, args.out_filename, args.output_dir)

def proteome_downloader(proteome_id, filename='input_proteome.fasta', output_dir=os.getcwd(), format_ = "fasta") -> None:
    """Downloads proteome from uniprot database into output multifasta file
    param: proteome_id: uniprot unique proteome id, not case-sensitive
    param: output_dir: output directory (default: current directory)
    param: format: uniprot API required format (default:fasta)
    param: filename: output proteome filename (default: input_proteome.fasta)
    """
    # use updated uniprot api:
    url=f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format={format_}&query=%28proteome%3A{proteome}%29'
    #url = f'https://www.uniprot.org/uniprot/?query=proteome:{proteome_id}&format={format}' 
    #filename = filename
    #output_dir = output_dir if output_dir != None else os.getcwd()
    response = requests.get(url, stream = True)
    text_file = open(os.path.join(output_dir, filename), 'wb')
    for chunk in response.iter_content(chunk_size=1024):
          text_file.write(chunk)
    # raise an AssertionError if the given proteome ID is not valid
    assert text_file.tell() > 0, f'{proteome_id} is not a valid Uniprot id'
    text_file.close()
    #print(f'Proteome {proteome_id} downloaded succesfully')
    #output = open(os.path.join(output_dir, filename), 'r').readlines()
    #return output
    return None
    
if __name__ == "__main__":
    main()
