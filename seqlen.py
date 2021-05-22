#!/usr/bin/env python

"""
Calculate sequence lengths in a multiFASTA file.
This script was inspired by Lesley Sitter's code published on biostars (www.biostars.org/p/148815/).

Command: python seqlen.py [input FASTA file] > [output TSV file]

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
"""

import sys

with open(sys.argv[1], "r") as f:
    input_fasta = f.read().splitlines()  # Newline characters are dropped.
print("\t".join(["Name", "Length"]))  # The header line
seq = ""
allow_write = False  # A flag indicating that the first sequence has been completely loaded.

for line in input_fasta:
    if line.startswith(">"):  # A new sequence is encountered
        if allow_write:
            print("\t".join([seqid, str(len(seq))]))
            seq = ""
        seqid = line.split(" ")[0]  # Get the sequence ID and ignore the sequence annotation
        seqid = seqid[1 : ]  # Drop the ">" character
        allow_write = True
    else:
        seq += line  # concatenate lines of the current sequence

# Write the length of the last sequence
print("\t".join([seqid, str(len(seq))]))
