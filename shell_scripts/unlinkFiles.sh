#!/bin/bash
# Iterates the unlink command on multiple soft links.
# Example:
#	unlinkFiles.sh reads/*.fastq.gz
#	unlinkFiles.sh $(cat link_list.txt)
# Author: Yu Wan (14-15 April 2017) published at https://github.com/wanyuac/BINF_toolkit
# License: Apache-2.0

links=( $@ )
for i in ${links[@]}; do
    unlink $i
done
