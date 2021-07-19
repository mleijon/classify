#!/usr/bin/python3

import argparse
import os
import gzip
from collections import defaultdict


class TxidHits:
    merged_results = list(dict(dict()))
    seq_set = set()

    def __init__(self, sa, txid):
        def make_reduced_daa(ddaf):
            red_daa = defaultdict(dict)
            seqname_set = set()
            red_seqs = dict()
            with gzip.open(ddaf, 'rt') as a:
                for line in a:
                    if line.split()[1] == txid:
                        red_daa[line.split(';')[0]]['size'] = \
                            line.split(';')[1][5:]
                        red_daa[line.split(';')[0]]['evalue'] = \
                            line.split()[2].strip()
                        seqname_set.add(line.split(';')[0])
            with gzip.open(ddaf.replace('.daa', '_uq.fa'), 'rt') as f:
                x = 0
                while True:
                    line = f.readline()
                    if not line:
                        break
                    if line.split(';')[0][1:] in seqname_set:
                        x= x + int(red_daa[line.split(';')[0][1:]]['size'])
                        red_daa[line.split(';')[0][1:]]['seq'] = next(f)
            print(x)
        self.sample = sa
        self.txid = id
        dirs = []
        self.home = os.path.join(input_dir, self.sample)
        for split in os.listdir(self.home):
            dirs.append(os.path.join(self.home, split))
        for subdir in dirs:
            print(subdir)
            for file in os.listdir(subdir):
                if file.endswith(".daa.gz"):
                    make_reduced_daa(os.path.join(subdir, file))


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('inp', type=str, nargs='*', help='result dict + taxids')
    ARGS = PARSER.parse_args()
    # Check arguments
    input_dir = os.path.normpath(ARGS.inp[0])
    if not os.path.isdir(input_dir):
        exit('First argument should be the directory "classification" '
             'containing the classified sequence data')
    for item in ARGS.inp[1:]:
        try:
            int(item)
        except ValueError:
            exit('The taxonomic ids should all be integers')
    samples = list()
    for filename in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, filename)):
            samples.append(filename)

    # Usearch produce fasta-files with max 80 columns. Here the new-lines
    # separating sequence is removed to produce a fasta-file with sequences on a
    # single line

    for sample in samples:
        for split_dir in os.listdir(os.path.join(ARGS.inp[0], sample)):
            parent_dir = os.path.join(ARGS.inp[0], sample)
            for file in os.listdir(os.path.join(parent_dir, split_dir)):
                if file.endswith('.fasta.gz'):
                    print('Reformatting fasta-file: {}'.format(file))
                    parent_dir = os.path.join(parent_dir, split_dir)
                    with gzip.open(os.path.join(parent_dir, file), 'rt')\
                            as f_in, gzip.open(os.path.join(parent_dir,
                                file.replace('fasta', 'fa')), 'wt') as f_out:
                        new_line = ''
                        for line in f_in:
                            if line[0] == '>':
                                f_out.write(new_line + line)
                                new_line = '\n'
                            else:
                                f_out.write(line.strip())
                    os.remove(os.path.join(parent_dir, file))

    test = TxidHits(samples[0], ARGS.inp[1])