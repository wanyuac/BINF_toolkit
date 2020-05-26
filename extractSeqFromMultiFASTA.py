#!/usr/bin/env python

"""
This script extracts sequences from a multi-FASTA file using sequence IDs, which are defined in sequence headers.

Input: a multi-FASTA file from stdin. It can be a gene-feature file (.ffn) downloaded from the NCBI nucleotide database.
Argument: a comma-delimited string of target sequence IDs.

Usage:
    cat input.fna | python extractSeqFromMultiFASTA.py "gene1,gene2,...,geneN" > output.fna
    Or,
    targets=$(cat genes.txt)  # genes.txt contains a single comma-delimited line.
    cat input.fna | python extractSeqFromMultiFASTA.py $targets > output.fna

Author: Yu Wan (wanyuac@126.com, https://github.com/wanyuac)
Python version 2 and 3 compatible
License: GNU GPL 2.1
First edition: 5 July 2016; the latest edition: 25 May 2020
Previous name: extract_fasta_loci.py
"""

from __future__ import print_function
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    # read the list of locus tags
    try:
        """
        First, drop the newline character in any combinations of \r and \n.
        Otherwise, the last ID does not match to any sequence ID.
        """
        loci = sys.argv[1].rstrip("\r\n")
        loci = loci.split(",")  # Parse the string for target sequence IDs
    except ValueError:
        print("Error: missing argument. A comma-delimited string of sequence IDs is required.")
    
    for seq in SeqIO.parse(sys.stdin, "fasta"):  # read the input FASTA file from stdin
        if seq.id in loci:
            print(seq.format("fasta"))  # write the current sequence to the stdout

if __name__ == "__main__":
    main()
    