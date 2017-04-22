#!/bin/sh
# Makes symbolic links of any kind of files in accordance with a list of sample names as inputs from
# either a single-column text file or stdin.
# Format of the command line:
#   bash mkSymbolicLins.sh [suffix] [source directory] [output directory] sample_names.txt
# Notice: 1. source and output directories must not be the same; 2. no forward slash ("/") should be attached to directory names.
# Examples:
#   sh mkVCFLinks.sh '_snps.vcf' ~/data ~/links strain_names.txt
#   sh mkVCFLinks.sh '_snps.vcf' . ~/links strain_names.txt
#   cat strain_names.txt | sh mkVCFLinks.sh '_snps.vcf' ~/data ~/links
#   In all examples, symbolic links [strain name]__snps.vcf will be created under the directory ~/links.
# Limitation: every pair of original file and its symobolic link shares the same filename suffix. Hence users must separate them
# with different directories.
# Author: Yu Wan (20, 22 Apr 2017)
# Licence: Apache-2.0

while IFS= read -r id; do
  #ln -s ${2}/${id}${1} ${3}/${id}${1}
  printf ${2}/${id}${1} ${3}/${id}${1}  # for debugging
done < "${4:-/dev/stdin}"  # takes $4 if defined otherwise takes the stdin
