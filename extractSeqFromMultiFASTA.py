#!/usr/bin/env python

"""
This script extracts sequences from a multi-FASTA file using sequence IDs, which are defined in sequence headers. It
basically filters contigs from a multi-FASTA file based on sequence IDs. The script ignores sequence annotation that
is separated from the sequence ID by a white space.

Input: a multi-FASTA file from stdin. It can be a gene-feature file (.ffn) downloaded from the NCBI nucleotide database
or an assembly file comprised of several contig sequences.

Argument: a comma-delimited string of target sequence IDs.

Usage:
    cat input.fna | python extractSeqFromMultiFASTA.py "gene1,gene2,...,geneN" > output.fna
    cat input.fna | python extractSeqFromMultiFASTA.py "contig1,contig2,...,contigM" > output.fna
    Or,
    targets=$(cat seqIDs.txt)  # seqIDs.txt contains a single comma-delimited line.
    cat input.fna | python extractSeqFromMultiFASTA.py $targets > output.fna
    For SPAdes assemblies:
    targets=$(cat seqIDs.txt)
    cat input__scaffolds.fna | sed 's/_length_/ /g' | python extractSeqFromMultiFASTA.py $targets > input__scaffolds_subset.fna

Author: Yu Wan (wanyuac@126.com, https://github.com/wanyuac)
Python version 2 and 3 compatible
License: GNU GPL 2.1
First edition: 5 July 2016; the latest edition: 13 Sep 2021
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
    