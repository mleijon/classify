#!/usr/bin/python3

import argparse
import os
import gzip
from collections import defaultdict


class TxidHits:
    seq_set = set()

    def __init__(self, sa, txid):
        self.merged_results = defaultdict(dict)

        def make_reduced_seqs(ddaf):
            red_seqs = defaultdict(dict)
            seqname_set = set()
            with gzip.open(ddaf, 'rt') as a:
                for line in a:
                    if line.split()[1] == txid:
                        red_seqs[line.split(';')[0]]['size'] = \
                            line.split(';')[1][5:]
                        red_seqs[line.split(';')[0]]['evalue'] = \
                            line.split()[2].strip()
                        seqname_set.add(line.split(';')[0])
            with gzip.open(ddaf.replace('.daa', '_uq.fa'), 'rt') as f:
                while True:
                    line = f.readline()
                    if not line:
                        break
                    if line.split(';')[0][1:] in seqname_set:
                        red_seqs[next(f).strip()] = \
                            red_seqs.pop(line.split(';')[0][1:])
            return red_seqs
        self.sample = sa
        self.txid = id
        dirs = []
        self.home = os.path.join(input_dir, self.sample)
        for split in os.listdir(self.home):
            dirs.append(os.path.join(self.home, split))
        for subdir in dirs:
            for file in os.listdir(subdir):
                if file.endswith(".daa.gz"):
                    split = make_reduced_seqs(os.path.join(subdir, file))
                    for key in split:
                        if key not in self.merged_results.keys():
                            self.merged_results.update({key:split[key]})
                        else:
                            new_size = str(int(self.merged_results[key]['size']
                                               ) + int(split[key]['size']))
                            self.merged_results[key]['size'] = new_size


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
    x = 0
    for key in test.merged_results:
        x += int(test.merged_results[key]['size'])
    print(x)
