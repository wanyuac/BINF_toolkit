#!/usr/bin/env python3
"""
Exclude nucleotide/protein sequences of pseudo genes from a multi-Fasta file.

Copyright (C) 2025 Yu Wan <wanyuac@gmail.com>
First release: 2 Jan 2025; latest update: 3 Jan 2025
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
"""

from Bio import SeqIO
from argparse import ArgumentParser
import re

def parse_arguments():
    parser = ArgumentParser(description = "Filter pseudo sequences from a multi-FASTA file.")
    parser.add_argument('--input', '-i', dest = 'input', required = True, help = "Path to the input FASTA file.")
    parser.add_argument('--output', '-o', dest = 'output', required = False, default = 'filtered_output.fasta', help = "Path to the output FASTA file.")
    parser.add_argument('--pseudo', '-p', dest = 'pseudo', required = False, default = 'pseudo.fasta', help = "Path to the output FASTA file of excluded sequences")
    parser.add_argument('--discard_annot', '-d', dest = 'discard_annot', required= False, action = 'store_true', help = "A flag to discard sequence annotations and only keep names")
    return parser.parse_args()

def filter_pseudo_sequences(input_fasta, output_fasta, pseudo_fasta, discard_annot):
    with open(input_fasta, 'r') as infile,\
        open(output_fasta, "w") as outfile,\
        open(pseudo_fasta, 'w') as pseudofile:
        for record in SeqIO.parse(infile, 'fasta'):  # Iterate through sequences in the input FASTA file
            if "[pseudo=true]" in record.description:  # Check if "[pseudo=true]" is in the header
                output_handle = pseudofile
            else:
                output_handle = outfile
            if discard_annot:
                record.id = rename_seq(record.description, record.id)
                record.description = record.id
            SeqIO.write(record, output_handle, 'fasta')

def rename_seq(seq_description, seq_id):
    match_locus_tag = re.search(r'\[locus_tag=([^\]]+)\]', seq_description)  # Extract locus_tag and protein_id from the description using regular expressions
    match_protein_id = re.search(r'\[protein_id=([^\]]+)\]', seq_description)
    locus_tag = match_locus_tag.group(1) if match_locus_tag else None
    protein_id = match_protein_id.group(1) if match_protein_id else None
    if locus_tag and protein_id:
        new_id = f"{locus_tag}__{protein_id}"
    else:
        new_id = seq_id  # No change to the sequence ID
    return new_id

def main():
    args = parse_arguments()
    filter_pseudo_sequences(args.input, args.output, args.pseudo, args.discard_annot)
    print(f"Filtered sequences have been written to {args.output} and {args.pseudo}.")

if __name__ == '__main__':
    main()
