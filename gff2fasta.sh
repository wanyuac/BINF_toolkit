#!/bin/bash
# This script extract sequence regions from GFF3 files which contain complete assembled sequences.
# In this kind of file, the sequences must be put at the end of each file. The sequence domain is
# separated from the annotation domain by the delimiter "###FASTA".
# Usage: bash gff2fasta.sh [input GFF file(s)]
# Examples:
#   bash gff2fasta.sh *.gff
#   bash gff2fasta.sh strain1.gff strain2.gff ...
# Licence: GNU GPL 2.1
# Author: Yu Wan (wanyuac@gmail.com)
# Development history: 21/7/2016

ext='fna'  # the file extension

for f in "$@"; do  # loop through each argument
    base=`basename $f .gff`  # remove the path as well as the file extension
    k=`grep -n '##FASTA' $f | cut -f1 -d ':'`
    tail -n +$((k + 1)) $f > ${base}.$ext  # print lines starting with the kth
done
