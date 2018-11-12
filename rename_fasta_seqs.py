#!/usr/bin/env python

"""
Rename sequences in a FASTA file. It filters out sequences that are not included in the target list,
when specified.

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Python version 2 and 3 compatible
License: GNU GPL 2.1
First edition: 11 Nov 2018, the latest revision: 13 Nov 2018.
Created and finished in Nara, Japan.
"""

from __future__ import print_function
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description="Read options and arguments")
    parser.add_argument("--fasta", "-f", dest = "fasta", type = str, required = True, help = "A FASTA file whose sequences will be renamed.")
    parser.add_argument("--mapping", "-m", dest = "mapping", type = str, required = True, help = "A tab-delimited file mapping original sequence IDs to new IDs.")
    parser.add_argument("--out", "-o", dest = "out", type = str, required = False, default = "./renamed.fasta", help = "Name and path for output FASTA file.")
    parser.add_argument("--keep_all", "-k", dest = "keep_all", action = "store_true", required = False, help = "Set to keep all sequences when some IDs are not found in the rename table.")

    return parser.parse_args()


def main():
    args = parse_arguments()
    mapping = import_mapping_table(args.mapping)
    to_rename = list(mapping.keys())
    in_fasta = open(args.fasta, "rU")
    out = open(args.out, "w")
    
    for seq in SeqIO.parse(in_fasta, "fasta"):  # read the input FASTA file from stdin
        if seq.id in to_rename:
            seq.id = mapping[seq.id]
            print(seq.format("fasta"), file = out)
        elif args.keep_all:
            print(seq.format("fasta"), file = out)

    in_fasta.close()
    out.close()
    
    return


def import_mapping_table(rename):
    # Read the tab-delimited table for renaming sequences.
    with open(rename, "rU") as f:
        lines = f.read().splitlines()
        
    r = {}
    for l in lines:
        old_id, new_id = l.split("\t")
        r[old_id] = new_id
        
    return(r)


if __name__ == "__main__":
    main()
