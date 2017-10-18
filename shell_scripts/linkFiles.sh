#!/bin/bash
# Creating symbolic links according to a table of two columns: original file path and link path, separated by commas.
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First and the latest edition: 18 Oct 2017

# Display help information ###############
display_usage(){
    echo "
    Usage:
           chmod a+x linkFiles.sh  # before the first run
           ./linkFiles.sh [input CSV file]
           cat [input CSV file] | ./linkFiles.sh
    The CSV file should not contain a header line. The first column consists of original file paths, and
    the second column consists of link paths.
    An example of the CSV file:
    ~/data/genome1_1.fasta,/scratch/input/genome1_unimelb.fna
    ~/data/genome1_2.fasta,/scratch/input/genome1_zju.fna
    
    Notice a user must ensure the directory is accessible for storing links.
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Otherwise, make symbolic links following the input file ###############
while read line; do
    IFS="," read -ra paths <<< "$line"  # split the delimited string into an arrary of two elements
    echo "ln -s ${paths[0]} ${paths[1]}"
done < "${1:-/dev/stdin}"  # accept either a file name or lines from stdin
