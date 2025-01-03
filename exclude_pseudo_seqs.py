#!/usr/bin/env python3
"""
Exclude nucleotide/protein sequences of pseudo genes from a multi-Fasta file.

Copyright (C) 2025 Yu Wan <wanyuac@gmail.com>
Release: 2 Jan 2025
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
"""

from Bio import SeqIO
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = "Filter pseudo sequences from a multi-FASTA file.")
    parser.add_argument('--input', '-i', dest = 'input', required = True, help = "Path to the input FASTA file.")
    parser.add_argument('--output', '-o', dest = 'output', required = False, default = 'filtered_output.fasta', help = "Path to the output FASTA file.")
    parser.add_argument('--pseudo', '-p', dest = 'pseudo', required = False, default = 'pseudo.fasta', help = "Patho to the output FASTA file of excluded sequences")
    return parser.parse_args()

def filter_pseudo_sequences(input_fasta, output_fasta, pseudo_fasta):
    with open(input_fasta, 'r') as infile,\
        open(output_fasta, "w") as outfile,\
        open(pseudo_fasta, 'w') as pseudofile:
        for record in SeqIO.parse(infile, 'fasta'):  # Iterate through sequences in the input FASTA file
            if "[pseudo=true]" in record.description:  # Check if "[pseudo=true]" is in the header
                SeqIO.write(record, pseudofile, 'fasta')
            else:
                SeqIO.write(record, outfile, 'fasta')

def main():
    args = parse_arguments()
    filter_pseudo_sequences(args.input, args.output, args.pseudo)
    print(f"Filtered sequences have been written to {args.output} and {args.pseudo}.")

if __name__ == '__main__':
    main()
