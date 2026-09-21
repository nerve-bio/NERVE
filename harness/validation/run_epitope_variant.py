#!/usr/bin/env python3
import argparse
import os
import pickle

from code.Epitope import epitope


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()

    with open(args.snapshot, 'rb') as handle:
        proteins = pickle.load(handle)

    os.makedirs(args.output)
    proteins = epitope(proteins, args.output, 9, 11, 1, 1, 0.9, 'False', 3,
                       args.workers)
    with open(os.path.join(args.output, '12_epitope.pkl'), 'wb') as handle:
        pickle.dump(proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
