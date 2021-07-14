#!/usr/bin/python3

from collections import Counter, defaultdict
import argparse
from pathlib import Path
home = str(Path.home())

ROOT = '1'
CELLULAR_ORGANISM = '131567'
ARCHAEA = '2157'
BACTERIA = '2'
EUCARYOTA = '2759'
VIRUSES = '10239'
UNCLASSIFIED = '12908'
OTHERSEQ = '28384'
tax_classes = [ARCHAEA, BACTERIA, EUCARYOTA, VIRUSES,
               UNCLASSIFIED, OTHERSEQ]
taxsets = {CELLULAR_ORGANISM: set(), ARCHAEA: set(), BACTERIA: set(),
           EUCARYOTA: set(), VIRUSES: set(), UNCLASSIFIED: set(),
           OTHERSEQ: set()}
taxcounts = defaultdict(int)
results = defaultdict(lambda: [])

# THese files are needed and downloaded from ncbi:
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

names_file = home + '/classify/names.dmp'
# Taxonomy names file has these fields:
#
# tax_id					-- the id of node associated with this name
# name_txt		     		-- name itself
# unique name				-- the unique variant of this name if name not
#                              unique
# name class				-- (synonym, common name, ...)

deleted_file = home + '/classify/delnodes.dmp'
# Deleted nodes (nodes that existed but were deleted) file field:
#
# tax_id					-- deleted node id

merged_file = home + '/classify/merged.dmp'
# Merged nodes file fields:
#
# old_tax_id                -- id of nodes which has been merged
# new_tax_id                -- id of nodes which is result of merging


lineage_file = home + '/classify/taxidlineage.dmp'


# Taxonomy id lineage file fields:
#
# tax_id                    -- node id
# lineage                   -- sequence of node ids separated by space denoting
#                              nodes' ancestors starting from the most distant
#                              one and ending with the immediate one


# Count the sequence reads (rows) in the input diamond file (.daa)
def count_rows(fname):
    with open(fname) as f:
        nrofreads = 0
        for i, l in enumerate(f):
            if ARGS.derep:
                nrofreads += int(l.split(';size=')[1].split(';\t')[0])
    if ARGS.derep:
        i = nrofreads
        return i
    return i + 1


# Counts the number of reads without hits, '0' in .daa file, and the number of
# reads assigned to 'ROOT', '1' in .daa file. The interpretation of '1' as
# a taxonomy id is that the reads have hit to the database but of equal
# similarity to enrries in the database of widely different taxa.
def count_0and1(fname):
    count0 = 0
    count1 = 0
    with open(fname) as f:
        for row in f:
            if ARGS.derep:
                nrofreads = int(row.split(';size=')[1].split(';\t')[0])
            else:
                nrofreads = 1
            item = row.split()
            if item[1].strip() == '0':
                count0 += nrofreads
            elif item[1].strip() == '1':
                count1 += nrofreads
    return count0, count1


# Create sets with the taxonomic classes: ARCHAEA, BACTERIA, EUCARYOTA, VIRUSES,
# UNCLASSIFIED, OTHERSEQ
def create_taxsets(fname):
    with open(fname) as f:
        for line in f:
            line = line.replace('|', ' ')
            item = list(map(str.strip, line.split()[0:3]))
            for taxclass in tax_classes:
                if taxclass in item:
                    taxsets[taxclass].add(item[0])
                    break
    return


def count_txclass(fname):
    with open(fname) as f:
        for row in f:
            taxid = row.split()[1].strip()
            if ARGS.derep:
                nrofreads = int(row.split(';size=')[1].split(';\t')[0])
            else:
                nrofreads = 1
            if taxid in taxsets[ARCHAEA]:
                taxcounts[ARCHAEA] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[ARCHAEA]):
                taxcounts[ARCHAEA] += nrofreads
            elif taxid in taxsets[BACTERIA]:
                taxcounts[BACTERIA] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[BACTERIA]):
                taxcounts[BACTERIA] += nrofreads
            elif taxid in taxsets[EUCARYOTA]:
                taxcounts[EUCARYOTA] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[EUCARYOTA]):
                taxcounts[EUCARYOTA] += nrofreads
            elif taxid in taxsets[VIRUSES]:
                taxcounts[VIRUSES] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[VIRUSES]):
                taxcounts[VIRUSES] += nrofreads
            elif taxid in taxsets[UNCLASSIFIED]:
                taxcounts[UNCLASSIFIED] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[UNCLASSIFIED]):
                taxcounts[UNCLASSIFIED] += nrofreads
            elif taxid in taxsets[OTHERSEQ]:
                taxcounts[OTHERSEQ] += nrofreads
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[OTHERSEQ]):
                taxcounts[OTHERSEQ] += nrofreads
            elif taxid in deletedTaxids:
                taxcounts['DELETED'] += nrofreads
            elif taxid == CELLULAR_ORGANISM:
                taxcounts[CELLULAR_ORGANISM] += nrofreads
    return


def parse_files(m, d, n):
    # Creates dictionary with merged taxids
    with open(m) as f:
        for line in f:
            item = line.split('|')
            mergedTaxids.update({item[0].strip(): item[1].strip()})

    # Creates a set with deleted taxids
    with open(d) as f:
        for line in f:
            item = line.split('|')[0].strip()
            deletedTaxids.add(item)

    # Creates a dictionary with key:values equal taxid:scientific name
    with open(n) as f:
        for line in f:
            item = line.split('|')
            if "scientific" in item[3]:
                sciNames.update({item[0].strip(): item[1].strip()})


def write_summary(fname):
    with open(fname, 'w') as f:
        f.write('Total nr of reads: {:>31}\n'.format(nr_of_reads))
        f.write('\nNo hits to db: {:>35} ({} %)\n'.format(nr_of_nohits[0], (
            round(100 * (nr_of_nohits[0] / nr_of_reads), 2))))
        f.write('Could not be classified by db: {:>19} ({} %)\n'.format(
            nr_of_nohits[1], (round(100 * (nr_of_nohits[1] / nr_of_reads), 2))))
        f.write('\nClassified by db: {:>32} ({} %)\n'.format(
            sum(taxcounts.values()),
            (round(100 * (sum(taxcounts.values()) / nr_of_reads), 2))))
        f.write('------------------------------'
                '------------------------------\n')
        f.write('Taxid deleted from db: {:>27} ({} %)\n'.format(
            taxcounts['DELETED'], (
                round(100 * (taxcounts['DELETED'] / nr_of_reads), 2))))
        f.write('CELLULAR ORGANISM: {:>31} ({} %)\n'.
                format(taxcounts[CELLULAR_ORGANISM],
                       (round(100 * (taxcounts[CELLULAR_ORGANISM]
                                     / nr_of_reads), 2))))
        f.write('ARCHAEA: {:>41} ({} %)\n'.format(taxcounts[ARCHAEA], (
            round(100 * (taxcounts[ARCHAEA] / nr_of_reads), 2))))
        f.write('BACTERIA: {:>40} ({} %)\n'.format(taxcounts[BACTERIA], (
            round(100 * (taxcounts[BACTERIA] / nr_of_reads), 2))))
        f.write('EUCARYOTA: {:>39} ({} %)\n'.format(taxcounts[EUCARYOTA], (
            round(100 * (taxcounts[EUCARYOTA] / nr_of_reads), 2))))
        f.write('VIRUS: {:>43} ({} %)\n'.format(taxcounts[VIRUSES], (
            round(100 * (taxcounts[VIRUSES] / nr_of_reads), 2))))
        f.write('UNCLASSIFIED: {:>36} ({} %)\n'.format(
            taxcounts[UNCLASSIFIED], (round(100 * (
                    taxcounts[UNCLASSIFIED] / nr_of_reads), 2))))
        f.write('OTHERSEQ: {:>40} ({} %)\n'.format(taxcounts[OTHERSEQ], (
            round(100 * (taxcounts[OTHERSEQ] / nr_of_reads), 2))))


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input diamond file',
                        required=True)
    PARSER.add_argument('-s', action='store_true',
                        help='switch for summary file output')
    PARSER.add_argument('-a', action='store_true',
                        help='switch for archaea file output')
    PARSER.add_argument('-b', action='store_true',
                        help='switch for bacteria file output')
    PARSER.add_argument('-e', action='store_true',
                        help='switch for eucaryota file output')
    PARSER.add_argument('-v', action='store_true',
                        help='switch for virus file output')
    PARSER.add_argument('-u', action='store_true',
                        help='switch for unclassified file output')
    PARSER.add_argument('-o', action='store_true',
                        help='switch for otherseq file output')
    PARSER.add_argument('-d', action='store_true',
                        help='switch for deleted file output')
    PARSER.add_argument('--derep', action='store_true',
                        help='switch for dereplicated reads')
    ARGS = PARSER.parse_args()
    sciNames = dict()
    mergedTaxids = dict()
    deletedTaxids = set()
    Results = []
    print('counting reads')
    nr_of_reads = count_rows(ARGS.f)
    print('counting \'no hits\' and \'root hits\'')
    nr_of_nohits = count_0and1(ARGS.f)
    print('loading taxonomy')
    create_taxsets(lineage_file)
    parse_files(merged_file, deleted_file, names_file)
    print('counting reads')
    count_txclass(ARGS.f)

    with open(ARGS.f) as f:
        for line in f:
            nrofreads = int(line.split(';size=')[1].split(';\t')[0])
            item = line.split()[1].strip()
            if item != '0':
                if item in mergedTaxids.keys():
                    if mergedTaxids[item] in taxsets[ARCHAEA]:
                        results[ARCHAEA] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                    elif mergedTaxids[item] in taxsets[BACTERIA]:
                        results[BACTERIA] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                    elif mergedTaxids[item] in taxsets[EUCARYOTA]:
                        results[EUCARYOTA] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                    elif mergedTaxids[item] in taxsets[VIRUSES]:
                        results[VIRUSES] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                    elif mergedTaxids[item] in taxsets[UNCLASSIFIED]:
                        results[UNCLASSIFIED] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                    elif mergedTaxids[item] in taxsets[OTHERSEQ]:
                        results[OTHERSEQ] += nrofreads*[
                            sciNames[mergedTaxids[item]]]
                elif item in deletedTaxids:
                    results['DELETED'] += [item + ' :deleted']
                else:
                    if item in taxsets[ARCHAEA]:
                        results[ARCHAEA] += nrofreads*[sciNames[item]]
                    elif item in taxsets[BACTERIA]:
                        results[BACTERIA] += nrofreads*[sciNames[item]]
                    elif item in taxsets[EUCARYOTA]:
                        results[EUCARYOTA] += nrofreads*[sciNames[item]]
                    elif item in taxsets[VIRUSES]:
                        results[VIRUSES] += nrofreads*[sciNames[item]]
                    elif item in taxsets[UNCLASSIFIED]:
                        results[UNCLASSIFIED] += nrofreads*[sciNames[item]]
                    elif item in taxsets[OTHERSEQ]:
                        results[OTHERSEQ] += nrofreads*[sciNames[item]]
    basename = ARGS.f.split('.daa')[0]
    if ARGS.a:
        with open(basename + '_archaea.txt', 'w') as f:
            f.write('{:<80}{}'.format('ARCHAEA', 'Nr of reads\n'))
            for item in Counter(results[ARCHAEA]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.b:
        with open(basename + '_bacteria.txt', 'w') as f:
            f.write('{:<80}{}'.format('BACTERIA', 'Nr of reads\n'))
            for item in Counter(results[BACTERIA]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.e:
        with open(basename + '_eucaryota.txt', 'w') as f:
            f.write('{:<80}{}'.format('EUCARYOTA', 'Nr of reads\n'))
            for item in Counter(results[EUCARYOTA]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.v:
        with open(basename + '_viruses.txt', 'w') as f:
            f.write('{:<80}{}'.format('VIRUSES', 'Nr of reads\n'))
            for item in Counter(results[VIRUSES]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.u:
        with open(basename + '_unclassified.txt', 'w') as f:
            f.write('{:<80}{}'.format('UNCLASSIFIED', 'Nr of reads\n'))
            for item in Counter(results[UNCLASSIFIED]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.o:
        with open(basename + '_otherseq.txt', 'w') as f:
            f.write('{:<80}{}'.format('OTHERSEQ', 'Nr of reads\n'))
            for item in Counter(results[OTHERSEQ]).most_common():
                f.write('{:<80}{}\n'.format(item[0], item[1]))
    if ARGS.d:
        with open(basename + '_deleted.txt', 'w') as f:
            f.write('{:<80}{}'.format('DELETED', 'Nr of reads\n'))
            for item in Counter(results['DELETED']).most_common():
                f.write('{:<80}{}\n'.format(item[0].split(':')[0].strip(),
                                            item[1]))
    if ARGS.s:
        summary_file_name = basename + '_summary.txt'
        write_summary(summary_file_name)
