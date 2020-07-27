#!/usr/bin/env python
"""
This script is derived from script extractSeqFromMultiFASTA.py and it aims to assign sequences
from an input multi-FASTA file into two multi-FASTA files based on a list of sequence IDs.

Usage:
    python splitMultiFASTA.py -i all_sequences.fna -l ids.txt -o1 in_list.fna -o2 out_list.fna

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 27 July 2020
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "Split a multi-FASTA file into two files")
    parser.add_argument("-i", type = str, required = True, help = "Input multi-FASTA file")
    parser.add_argument("-l", type = str, required = True, help = "Input list of sequence IDs for inclusion. One ID per line.")
    parser.add_argument("-o1", type = str, required = True, help = "Output FASTA file of sequences whose IDs comprise the ID list.")
    parser.add_argument("-o2", type = str, required = True, help = "Output FASTA file of sequences whose IDs are not on the ID list.")

    return parser.parse_args()


def main():
    args = parse_arguments()
    ids_inc = read_ID_list(args.l)
    
    out_1 = open(args.o1, "w")
    out_2 = open(args.o2, "w")
    for seq in SeqIO.parse(args.i, "fasta"):  # read the input FASTA file from stdin
        if seq.id in ids_inc:
            print(">" + seq.description, file = out_1)
            print(seq.seq, file = out_1)
        else:
            print(">" + seq.description, file = out_2)
            print(seq.seq, file = out_2)
    out_1.close()
    out_2.close()


def read_ID_list(list_file):
    with open(list_file, "r") as f:
        ids = f.read().splitlines()

    return ids


if __name__ == "__main__":
    main()
    