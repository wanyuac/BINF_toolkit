#!/bin/bash
# This script concatenates reference sequences of the same bacterial strain into a single multi-FASTA file.
# The name of every FASTA file must follow the format: [strain name]__[accession number].fasta
# For example, AH0650_Sm1__LFJS01000001.fasta.
# Example command line:
#   bash concat_fasta.sh 'fasta/*__*.fasta'  # Quotes are necessary!
#   bash concat_fasta.sh 'fasta/strain__*.fasta'
#	bash concat_fasta.sh 'fasta/strain__*.faa'
# Licence: GNU GPL 2.1
# Author: Yu Wan (wanyuac@gmail.com)
# Development history: 9 Aug 2016, 12 Sep 2016

f=$1
ext=${f##*.}  # get the file name extension

strains=$(ls -1 ${1} | xargs -I '{}' basename {} ".${ext}" | grep -oP '.+(?=__)' | sort -u)
echo "$(echo $strains | tr " " "\n" | wc -l) strains are to be processed."

path=$(dirname "$1")
for s in ${strains}; do
    cat ${path}/${s}__*.${ext} > ${s}.${ext}
done
