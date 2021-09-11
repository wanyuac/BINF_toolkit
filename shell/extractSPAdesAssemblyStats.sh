#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 16 Apr 2021; latest update: 11 Sep 2021

display_usage() {
    echo "
    Extract contig/scaffold names, lengths, and depths from a SPAdes output FASTA file and save them in a tab-delimited text file.
    Command line:
      extractSPAdesAssemblyStats.sh input.fasta > output.tsv  # Default: do not print a header line (for the convenience of concatenating files)
      extractSPAdesAssemblyStats.sh input.fasta 1 > output.tsv  # Print a header line"
}

if [ -z $1 ]; then
    display_usage
    exit
fi

if [ "$2" -eq "1" ]; then
    echo -e 'Node\tLength\tDepth'  # Print the header line. Note that the echo command automatically appends a newline character to the output string.
fi

grep '>' $1 | sed -e 's/>//g' | sed -e 's/_length_/\t/g' | sed -e 's/_cov_/\t/g'