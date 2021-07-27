#!/usr/bin/python3

import argparse
import os
import re
from fasta import FastaList


def create_filter(seq):
    seq = seq.lower()
    filt = ''
    for char in seq:
        if char in ['a', 'c', 'g', 't']:
            filt += char
        if char == 'y':
            filt += '[ct]'
        if char == 'r':
            filt += '[ag]'
        if char == 'k':
            filt += '[gt]'
        if char == 'm':
            filt += '[ac]'
        if char == 's':
            filt += '[cg]'
        if char == 'w':
            filt += '[at]'
        if char == 'b':
            filt += '[cgt]'
        if char == 'd':
            filt += '[agt]'
        if char == 'h':
            filt += '[act]'
        if char == 'v':
            filt += '[acg]'
        if char in ['x', 'n']:
            filt += '[acgt]'
    return filt


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='New filter fasta-file will be'
                                                 ' created containing only '
                                                 'sequence with the sub-'
                                                 'sequences listed in a filter-'
                                                 'file. The directory given by'
                                                 'the -d option should contain'
                                                 'fasta-files from '
                                                 'taxid2seqs.py and a filter '
                                                 'file named "filter". The '
                                                 'output files wil be tagged'
                                                 'with "_filt"')
    PARSER.add_argument('-d', type=str, help='directory containing input data'
                                             'and teh "filter"-file',
                        required=True)
    ARGS = PARSER.parse_args()
    # Check files
    if os.path.isfile(os.path.join(ARGS.d, 'filter.fa')):
        primers = FastaList(os.path.join(ARGS.d, 'filter.fa'))
    else:
        exit('filter-file not found')
    samples = list()
    for file in os.listdir(ARGS.d):
        if file.endswith("sel_txids.fasta"):
            samples.append(file.split('_')[0])
    if not samples:
        exit('No input fasta-files found')
    primer_set1 = [primers.seq_list[0].split()[1].lower()]
    primer_set1 += [primers.seq_list_revc()[1].split()[1].lower()]
    patterns1 = []
    primer_set2 = [primers.seq_list[1].split()[1].lower()]
    primer_set2 += [primers.seq_list_revc()[0].split()[1].lower()]
    patterns2 = []
    for primer in primer_set1:
        patterns1 += [re.compile(create_filter(primer))]
    print(patterns1)

