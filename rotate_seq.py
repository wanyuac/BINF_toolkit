#!/usr/bin/env python
"""
This script rotates circular sequences and makes them to start from specific positions in output FASTA files. This script is
useful for improving complete genome assemblies as well as read simulation.

Command:
    python rotate_seq.py -i input.fna -t new_starts.tsv 1> output.fna 2> messages.err
    python rotate_seq.py -i input.fna -t new_starts.tsv | gzip > output.fna.gz  # Only stdout goest into the pipe.

New start positions are specified in a two-column, tab-delimited, header-free table (parameter '-t'): sequence ID\tPosition.
No change is applied to input sequences whose sequence IDs are not listed in this table (e.g., when some sequences are linear
or incomplete).

Dependencies: Python 3, BioPython 1.78+

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 17 June 2021; the latest update: 18 June 2021
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq  # Bio.Alphabet has been removed from BioPython from v1.78. See https://biopython.org/wiki/Alphabet.
from argparse import ArgumentParser


def parse_argument():
    parser = ArgumentParser(description = "Restart circular sequences from given positions")
    parser.add_argument("-i", "--input", dest = "i", type = str, required = True, help = "An input FASTA file")
    parser.add_argument("-t", "--table", dest = "t", type = str, required = True, \
    help = "A tab-delimited, header-free table of two columns: sequence ID, the base to be used as the first base of the new sequence")
    return parser.parse_args()


def main():
    args = parse_argument()
    prev_seqs = import_seqs(args.i)
    pos_spec = import_positions(args.t)
    for i, contig in prev_seqs.items():
        if i in pos_spec.keys():
            p = pos_spec[i]
            if p > 0:  # Convert the p-th position into Python's character index. Note that no change will be carried out by Python if p > len(contig.seq).
                ori = p - 1
            else:
                print("Warning: position %i for sequence %s cannot be negative. No change will be applied to this sequence." % (p, i), file = sys.stderr)
                ori = 0
            print("Restart sequence %s from base %i" % (i, p), file = sys.stderr)
            s = str(contig.seq)
            contig.seq = Seq(s[ori : ] + s[ : ori])  # Rotate the current sequence; no change when p = 0. "generic_dna" is no longer needed from BioPython v1.78.
        else:
            print("Warning: sequence " + i + " is not found in the position table. No change will be applied to this sequence.", file = sys.stderr)
        SeqIO.write(contig, sys.stdout, "fasta")
    return


def import_seqs(fasta):
    check_file(fasta)
    seqs = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    return seqs


def import_positions(tsv):
    check_file(tsv)
    with open(tsv, "r") as f:
        tab = dict()  # {id : new start position}
        lines = f.read().splitlines()
        for line in lines:
            i, p = line.split("\t")
            tab[i] = int(p)
    return(tab)


def check_file(f):
    if not os.path.exists(f):
        print("Argument error: file " + f + " is not accessible.", file = sys.stderr)
        sys.exit(1)
    return


if __name__ == "__main__":
    main()
