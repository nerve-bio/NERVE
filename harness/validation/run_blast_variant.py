#!/usr/bin/env python3
import argparse
import os
import pickle

from code.Immunity import autoimmunity, conservation, mouse


def save(path, proteins):
    with open(path, 'wb') as handle:
        pickle.dump(proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)


parser = argparse.ArgumentParser()
parser.add_argument('--snapshot', required=True)
parser.add_argument('--proteome1', required=True)
parser.add_argument('--proteome2', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--threads', type=int, required=True)
args = parser.parse_args()

with open(args.snapshot, 'rb') as handle:
    proteins = pickle.load(handle)

os.makedirs(args.output)
proteins = autoimmunity(proteins, args.proteome1, args.output,
                        '/usr/nerve_python/NERVE', 1e-10, 9, 1, 3,
                        threads=args.threads)
save(os.path.join(args.output, '06_autoimmunity.pkl'), proteins)
proteins = mouse(proteins, args.output, '/usr/nerve_python/NERVE', 1e-10,
                 args.proteome1, 9, 3, 1, threads=args.threads)
save(os.path.join(args.output, '07_mouse.pkl'), proteins)
os.makedirs(os.path.join(args.output, 'compare_proteome'))
proteins = conservation(proteins, args.output, '/usr/nerve_python/NERVE',
                        1e-10, args.proteome1, args.proteome2, 9, 3, 1,
                        threads=args.threads)
save(os.path.join(args.output, '08_conservation.pkl'), proteins)
