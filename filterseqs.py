#!/usr/bin/python3

import argparse
import os
import re
from fasta import FastaList


def create_filter(seq, direction):
    min_length = ARGS.m
    opt_count = len(seq) - min_length
    seq = seq.lower()
    filt = ''
    count = 0
    for char in seq:
        count += 1
        if char in ['a', 'c', 'g', 't']:
            filt += char
        elif char == 'y':
            filt += '[ct]'
        elif char == 'r':
            filt += '[ag]'
        elif char == 'k':
            filt += '[gt]'
        elif char == 'm':
            filt += '[ac]'
        elif char == 's':
            filt += '[cg]'
        elif char == 'w':
            filt += '[at]'
        elif char == 'b':
            filt += '[cgt]'
        elif char == 'd':
            filt += '[agt]'
        elif char == 'h':
            filt += '[act]'
        elif char == 'v':
            filt += '[acg]'
        elif char in ['x', 'n']:
            filt += '[acgt]'
        else:
            exit('Error in primer sequence')
        if direction == 'fwd':
            if count <= opt_count:
                filt += '?'
            else:
                pass
        elif direction == 'rev':
            if count > min_length:
                filt += '?'
            else:
                pass
        else:
            exit('Primer direction error')
    return filt


def make_pattern_sets():
    primers = FastaList(os.path.join(ARGS.d, 'filter.fa'))
    primer_set1 = [primers.seq_list[0].split()[1].lower()]
    primer_set1 += [primers.seq_list_revc()[1].split()[1].lower()]
    patter_set1 = []
    direction = 'fwd'
    for primer in primer_set1:
        patter_set1 += [re.compile(create_filter(primer, direction))]
        direction = 'rev'
    primer_set2 = [primers.seq_list[1].split()[1].lower()]
    primer_set2 += [primers.seq_list_revc()[0].split()[1].lower()]
    pattern_set2 = []
    direction = 'fwd'
    for primer in primer_set2:
        pattern_set2 += [re.compile(create_filter(primer, direction))]
        direction = 'rev'
    pattern_sets = []
    pattern_sets += [patter_set1]
    pattern_sets += [pattern_set2]
    return pattern_sets


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
    PARSER.add_argument('-m', type=int, help='The minimum nr of matching'
                                             'nucleotides required for a match',
                        default=10, required=False)
    ARGS = PARSER.parse_args()
    # Check files
    if os.path.isfile(os.path.join(ARGS.d, 'filter.fa')):
        patterns = make_pattern_sets()
    else:
        exit('filter-file not found')
    samples = list()
    for file in os.listdir(ARGS.d):
        if file.endswith("sel_txids.fasta"):
            samples.append(file.split('_')[0])
    if not samples:
        exit('No input fasta-files found')
    else:
        for sample in samples:
            inp_file = os.path.join(ARGS.d, sample + '_sel_txids.fasta')
            for item in FastaList(inp_file).seq_list:
                seq = item.split()[1].lower()
                for pattern_set in patterns:
                    for primer_pattern in pattern_set:
                        print(primer_pattern)
                        print(seq)
                        print(primer_pattern.search(seq))
