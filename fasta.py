#!/usr/bin/python3
"""General purpose module to handle a fasta file object. Converts fastq to
fasta and unzip gzipped files"""
import sys
import os
import copy


class FastaListIterator:
    def __init__(self, fastalist):
        self._fastalist = fastalist
        self._index = 0

    def __next__(self):
        if self._index < self._fastalist.nr_seq:
            self._index += 1
            return self._fastalist.seq_list[self._index - 1]
        raise StopIteration


class FastaList:
    """From a fasta text file (fastq or fasta) creates a list with the id:s of
    the fasta sequences (id_list) and a list of the sequences with the newlines
    removed (seq_list), suitable for searching. The sequenes list elements are
    strings with a newline separating the id and the nucleotides sequence and
    ends with \n. id_list strings end without \n.
    """

    def __init__(self, fasta_name):  # Initialized by a filename string
        self.name = str(os.path.abspath(fasta_name))
        if self.name.upper().endswith('.GZ'):
            import gzip
            fa_file = open(self.name[:-3], 'w')
            with gzip.open(self.name, 'rt') as in_fi:
                fa_file.write(in_fi.read())
            fa_file.close()
            self.name = self.name[:-3]
        if self.name.upper().endswith(('.FQ', '.FASTQ')):
            self.fa_file = self.fq2fa()  # Handle file as fastq
        elif self.name.upper().endswith(('.FA', '.FASTA', 'FNA', 'AFA')):
            self.fa_file = self.rdfi()  # Handle file as fasta
        else:
            sys.exit('Not a fastq or fasta file')
        self.seq_list = []
        self.id_list = []
        newseq = ''
        for line in self.fa_file:
            if line.startswith('>'):
                if self.id_list:
                    self.seq_list.append(newseq + '\n')
                self.id_list.append(line[1:].strip())
                newseq = line
            else:
                newseq += line.strip()
        self.seq_list.append(newseq + '\n')
        self.fa_file.seek(0)
        self.nr_seq = len(self.id_list)

    def __iter__(self):
        return FastaListIterator(self)

    def rdfi(self):
        """ Return a fasta file object"""
        self.fa_file = open(self.name)
        return self.fa_file

    def fq2fa(self):
        """ Reads and return a fastq file converted to fasta format"""
        out_fa = open(self.name[:self.name.find('.')] + '.fa', 'w')
        for line in enumerate(self.rdfi()):
            if line[1][0] == '@' and line[0] % 4 == 0:
                out_fa.write(line[1].replace('@', '>', 1))  # seqid line
            elif line[0] % 4 == 1:
                out_fa.write(line[1])  # sequence line
            elif line[1][0] == '+' and line[0] % 4 == 2:
                pass
            elif line[0] % 4 == 3:
                pass
            else:
                out_fa.close()
                sys.exit('Not a fastq-file')
        out_fa.close()
        out_fa = open(out_fa.name)
        return out_fa

    def seq_list_rev(self):
        """Returns the reverse of a nucleotide sequence"""
        rev_seq_list = []
        for item in self.seq_list:
            seqid = item.split('\n')[0] + '_rev\n'
            seq = item.split('\n')[1][::-1] + '\n'
            rev_seq_list.append(seqid + seq)
        return rev_seq_list

    def seq_list_revc(self):
        """Creates a list the reverse complement of a nucleotide sequences in
        the list seq_list"""

        def comp(nt):
            return {
                'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G',
                'Y': 'R',
                'R': 'Y',
                'K': 'M',
                'M': 'K',
                'S': 'S',
                'W': 'W',
                'N': 'N',
                'B': 'V',
                'V': 'B',
                'D': 'H',
                'H': 'D',
                '-': '-'
            }[nt.upper()]

        seq_list_rc = []
        for item in self.seq_list:
            seqid = item.split('\n')[0] + '_rev\n'
            seq = ''
            for nt in item.split('\n')[1][::-1]:
                seq += comp(nt)
            seq_list_rc.append(seqid + seq)
        return seq_list_rc

    def divide(self, divisor):
        tmp_list = self.seq_list.copy()
        seq_list_div = []
        item_per_lst = self.nr_seq // divisor + (self.nr_seq % divisor)\
            // divisor
        for i in range(divisor - 1):
            seq_list_div.append(tmp_list[0:item_per_lst])
            del tmp_list[0:item_per_lst]
        seq_list_div.append(tmp_list)
        return seq_list_div

    def rmprimers(self, primer_file_name, primer_frac=0.5):
        # Removes primer sequences (listed in "primer_file_name) including
        # leading and trailing nucletides. Only the fraction primer_frac from
        # innermost part of the primer needs to be present

        seqs_noprimers = []  # List of sequences without primers

        # Create list of primers - primerlst
        pr = FastaList(primer_file_name)
        primers = pr.seq_list
        primers_rc = FastaList(primer_file_name).seq_list_revc()
        primers_all = primers + primers_rc
        primerlst = []
        for item in primers_all:
            primerseq = item.split('\n')[1]
            # Shorten the primers to keep the fraction - primer_frac
            primerlst.append(
                primerseq[round(len(primerseq) * (1 - primer_frac)):])
        for seq in self.seq_list:
            seqid = seq.split('\n')[0]
            seqseq = seq.split('\n')[1]
            found_primers = 0
            for primer in primerlst:
                if found_primers > 0:
                    seqs_noprimers.append(seqid + '\n' + seqseq + '\n')
                    break
                if primer in seqseq:
                    if seqseq.find(primer) < (
                            len(seqseq) - seqseq.find(primer)):
                        seqseq = seqseq[(seqseq.find(primer) + len(primer)):]
                    else:
                        seqseq = seqseq[:seqseq.find(primer)]
                    found_primers += 1
        pr.rdfi().close()
        return seqs_noprimers

    def crop_ends(self):
        # Remove aligned fasta columns from both ends with gaps until a column
        # without gaps are found

        def crop(alignment):
            hyphen = True
            while hyphen:
                hyphen = False
                for seq in alignment:
                    if seq.split('\n')[1][0] == '-':
                        hyphen = True
                        for i in range(len(alignment)):
                            seqid = alignment[i].split('\n')[0]
                            newseq = alignment[i].split('\n')[1][1:]
                            alignment[i] = seqid + '\n' + newseq + '\n'
                        break
            return alignment

        new_alignment = copy.deepcopy(self.seq_list)
        new_alignment = crop(new_alignment)
        for i in range(len(new_alignment)):
            seqid = new_alignment[i].split('\n')[0]
            seqseq = new_alignment[i].split('\n')[1]
            new_alignment[i] = seqid + '\n' + seqseq[::-1] + '\n'
        new_alignment = crop(new_alignment)
        for i in range(len(new_alignment)):
            seqid = new_alignment[i].split('\n')[0]
            seqseq = new_alignment[i].split('\n')[1]
            new_alignment[i] = seqid + '\n' + seqseq[::-1] + '\n'
        return new_alignment

    def wr_fasta_file(self, fasta_file_name, primer_file=None, mode='w'):
        fi = open(fasta_file_name, mode)
        if not primer_file:
            for item in self.seq_list:
                fi.write(item)
        else:
            for item in self.rmprimers(primer_file):
                fi.write(item)
        fi.close()

    def rm_non_agct_columns(self):
        indices2rm = set()
        allowed_nts = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', '-'}
        len_first_seq = len(self.seq_list[0].split('\n')[1])
        for item in self.seq_list:
            if len(item.split('\n')[1]) != len_first_seq:
                exit('sequence strings not of same length')
            for index, nt in enumerate(item.split('\n')[1]):
                if nt not in allowed_nts:
                    indices2rm.add(index)
        new_seq_list = []
        for item in self.seq_list:
            newseq = ''
            seqid, seq, _ = item.split('\n')
            for index, nt in enumerate(seq):
                if index in indices2rm:
                    continue
                else:
                    newseq += nt
            new_seq_list.append(seqid + '\n' + newseq + '\n')
            self.seq_list = new_seq_list


if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser(description='Convert fastq to fasta. '
                                                 'Assumes file extension'
                                                 ' "fastq"')
    PARSER.add_argument('-f', type=str, help='fastq filename', required=True)
    ARGS = PARSER.parse_args()
    fqfi = FastaList(ARGS.f)
    if ('.FQ' or '.FASTQ') in ARGS.f.upper():
        fqfi.fq2fa()
