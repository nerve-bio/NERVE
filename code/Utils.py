#!/usr/local/bin/python
"""NERVE useful functions and classes"""

import argparse, os, subprocess

def bashCmdMethod(bashCmd):
    """Run bash commands
    param: bashCmd: bash command to be run"""
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error
    
def dir_path(path:str) -> str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path
