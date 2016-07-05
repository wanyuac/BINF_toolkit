"""
This script extract loci from a FASTA file in accordance with a list of locus tags (namely, sequence names).
Number of options: one, which is a comma-delimited list of locus_tags in the GenBank file.
An error will arise if the list of locus tags is missing.

Usage:
    cat input.fna | python extract_fasta_loci.py "locus1,locus2,...,locusN" > output.fna
    Or,
    targets=$(cat locus_tags.txt)
    cat input.fna | python extract_fasta_loci.py $targets > output.fna

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Development history: 5 July 2016
Python version: 2.7.10
License: GNU GPL 2.1
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    # read the list of locus tags
    try:
        loci = sys.argv[1].split(",")
    except ValueError:
        print "Missing argument: a comma-separated list of locus tags is necessary."
    
    for seq in SeqIO.parse(sys.stdin, "fasta"):  # read the input FASTA file from stdin
        if seq.id in loci:
            print seq.format("fasta")  # write the current sequence to the stdout

if __name__ == "__main__":
    main()
    