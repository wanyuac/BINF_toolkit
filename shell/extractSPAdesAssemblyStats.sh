#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 16 Apr 2021; latest update: 11 Sep 2021

display_usage() {
    echo "
    Extract contig/scaffold names, lengths, and depths from a SPAdes output FASTA file and save them in a tab-delimited text file.
    Command line:
      extractSPAdesAssemblyStats.sh [input.fasta] > [isolate1.tsv]  # Single-assembly mode: print a header line
      extractSPAdesAssemblyStats.sh [input.fasta] [isolate name] > [isolate name.tsv]  # Multi-assembly mode: do not print a header line and append
      the isolate name in each line for the convenience of concatenating files. This mode is used in a loop that runs this script iteratively."
}

if [ -z $1 ]; then
    display_usage
    exit
fi

if [ -z "$2" ]; then  # Single-assembly mode
    echo -e 'Node\tLength\tDepth'  # Print the header line. Note that the echo command automatically appends a newline character to the output string.
    grep '>' $1 | sed -e 's/>//g' | sed -e 's/_length_/\t/g' | sed -e 's/_cov_/\t/g'
else  # Multi-assembly mode (namely, to loop through multiple FASTA files, where each iteration calls this script)
    IFS=$'\n'  # https://stackoverflow.com/questions/8768420/how-to-convert-command-output-to-an-array-line-by-line-in-bash
    lines=( $(grep '>' $1 | sed -e 's/>//g' | sed -e 's/_length_/\t/g' | sed -e 's/_cov_/\t/g') )
    for i in ${lines[@]}; do
        echo -e "${2}\t${i}"  # Add the assembly name to the head of the result line; use 'echo', not 'printf' (which does not print a newline character at the end of the output)
    done
fi
