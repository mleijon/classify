#!/usr/bin/python3

import argparse
import os
import regex as re
from fasta import FastaList


def create_filter(prm):
    prm = prm.lower()
    filt = ''
    count = 0
    for char in prm:
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
    filt = '(' + filt + ')'
    if ARGS.fuzzy_match:
        filt += '{s<=' + ARGS.fuzzy_match + '}'
    return filt


def make_pattern_sets():
    primers = FastaList(os.path.join(ARGS.d, 'filter.fa'))
    fwd_primer = primers.seq_list[0].split()[1]
    rev_primer = primers.seq_list[1].split()[1]
    if not ARGS.m:
        crop_fwd = 0
        crop_rev = 0
    else:
        pattern_size = ARGS.m
        if pattern_size >= min(len(fwd_primer), len(rev_primer)):
            exit('pattern-size larger than primer length')
        crop_fwd = len(fwd_primer) - pattern_size
        crop_rev = len(rev_primer) - pattern_size
    primer_set = [primers.seq_list[0].split()[1].lower()[crop_fwd:]]
    primer_set += [primers.seq_list_revc()[1].split()[1].lower()[crop_rev:]]
    pattern_pairs = []
    pattern_pair = []
    for primer in primer_set:
        pattern_pair += [re.compile(create_filter(primer))]
    pattern_pairs += [pattern_pair]
    if ARGS.reverse_patterns:
        pattern_pair = []
        primer_set = [primers.seq_list[1].split()[1].lower()[crop_rev:]]
        primer_set += [primers.seq_list_revc()[0].split()[1].lower()[crop_fwd:]]
        for primer in primer_set:
            pattern_pair += [re.compile(create_filter(primer))]
        pattern_pairs += [pattern_pair]
    return pattern_pairs


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
                                         'and the "filter"-file',
                    required=True)
PARSER.add_argument('-m', type=int, help='The minimum nr of matching'
                                         'nucleotides required for a match',
                    required=False)
# TBD
# PARSER.add_argument('--amplicon_fasta', action='store_true',
#                     help='switch for export of amplicon fasta files')
PARSER.add_argument('--hits', action='store_true',
                    help='switch for creating a primer hit list')
PARSER.add_argument('--reverse_patterns', action='store_true',
                    help='switch for including the reverse primer pattern')
PARSER.add_argument('--fuzzy_match', type=str,
                    help='switch for including the reverse primer pattern',
                    required=False)
ARGS = PARSER.parse_args()
# Check files
patterns = list()
if os.path.isfile(os.path.join(ARGS.d, 'filter.fa')):
    primer_list = FastaList(os.path.join(ARGS.d, 'filter.fa'))
    patterns = make_pattern_sets()
else:
    exit('filter-file not found')
samples = list()
for file in os.listdir(ARGS.d):
    if file.endswith("_filt.fasta"):
        os.remove(os.path.join(ARGS.d, file))
for file in os.listdir(ARGS.d):
    if file.endswith(".fasta"):
        samples.append(file.split('.')[0])
if not samples:
    exit('No input fasta-files found')
else:
    for sample in samples:
        inp_file = os.path.join(ARGS.d, sample + '.fasta')
        out_file = os.path.join(ARGS.d, sample + '_filt.fasta')

        if ARGS.hits:
            hit_file = os.path.join(ARGS.d, sample + '_hit.list')
            h = open(hit_file, 'w')
        filtered_seqs = set()
        with open(out_file, 'w') as f:
            for item in FastaList(inp_file).seq_list:
                seq = item.split('\n')[1].lower()
                seq_id = item.split('\n')[0]
                directions = set()
                for pattern_set in patterns:
                    hitlist = list()
                    direc = 'fwd'
                    for primer_pattern in pattern_set:
                        print(patterns)
                        print(pattern_set)
                        print(primer_pattern)
                        exit()
                        for p in primer_pattern.finditer(
                                seq, overlapped=False):
                            directions.add(direc)
                            hitlist.append((p.group(),
                                            p.fuzzy_counts[0],
                                            p.span()))
                        direc = 'rev'
                    if len(directions) == 2:
                        filtered_seqs.add(item)
                        if ARGS.hits:
                            h.write(seq_id + '\n')
                            for hit in hitlist:
                                h.write('{}: mismatches: {}; span: {}\n'
                                        .format(hit[0], hit[1], hit[2]))
            for item in filtered_seqs:
                f.write(item)
        if ARGS.hits:
            h.close()
