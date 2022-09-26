#!/usr/bin/env python

"""
Calculate sequence lengths in a multiFASTA file.
This script was inspired by Lesley Sitter's code published on biostars (www.biostars.org/p/148815/).

Command: python seqlen.py -i [input FASTA file] (-a) (-n) > [output TSV file]
Examples:
    python seqlen.py -i input.fna -a > seq_lengths.tsv  # With sequence annotation in the sequence description
    python seqlen.py -i input.fna -n > seq_lengths.tsv  # Only keep the sequence ID in the sequence description and ignore 'N' and '-' characters

Any character: a flag to keep sequence annotation in the output.

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Release: 2021; latest update: 26 Sep 2022
"""

import sys
import re
from argparse import ArgumentParser

def parse_argument():
    parser = ArgumentParser(description = "Calculating lengths of sequences in a FASTA file")
    parser.add_argument('-i', '--input', dest = 'i', type = str, required = True, help = "An input FASTA file")
    parser.add_argument('-a', '--annot', dest = 'a', action = 'store_true', help = "Keep sequence annotations in addition to sequence names")
    parser.add_argument('-n', '--nucl', dest = 'n', action = 'store_true', help = "Ignoring \'-\' and \'N\' in nucleotide sequences.")
    return parser.parse_args()

def main():
    args = parse_argument()

    with open(args.i, "r") as f:
        input_fasta = f.read().splitlines()  # Newline characters are dropped.
    print("\t".join(["Name", "Length"]))  # The header line
    seq = ""
    allow_write = False  # A flag indicating that the first sequence name has been completely loaded.
    ignore_annot = not args.a

    for line in input_fasta:
        if line.startswith(">"):  # A new sequence is encountered
            if allow_write:
                if args.n:
                    seq = re.sub('[N-]', '', seq.upper())  # Stripping multiple characters from a string. Ref: stackoverflow.com/questions/3900054/python-strip-multiple-characters.
                print("\t".join([seqid, str(len(seq))]))  # Write the name and length of the previous sequence
                seq = ""
            if ignore_annot:
                seqid = line.split(" ")[0]  # Get the sequence ID and ignore the sequence annotation
            else:
                seqid = line
            seqid = seqid[1 : ]  # Drop the ">" character
            allow_write = True
        else:
            seq += line  # concatenate lines of the current sequence

    print("\t".join([seqid, str(len(seq))]))  # Write the length of the last sequence
    return

if __name__ == '__main__':
    main()
