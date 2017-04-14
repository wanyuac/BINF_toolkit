#!/bin/bash
# Iterates the unlink command on multiple soft links.
# Example:
#   unlinkFiles.sh "*.fasta"
#	unlinkFiles.sh "reads/*.fastq.gz"
# Author: Yu Wan (14 April 2017) published at https://github.com/wanyuac/BINF_toolkit
# License: Apache-2.0

wildcard=$1
links=(`ls -1 ${wildcard}`)  # read a list of files into an array
for i in ${links[@]}; do
    unlink $i
done
