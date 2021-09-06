#!/bin/bash
# Creating symbolic links according to a table of two columns: original file path and link path, separated by tab characters.
# Copyright (C) 2017-2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First edition: 18 Oct 2017, the latest update: 6 Sep 2021
# Update note: changed the input file format from CSV to TSV for convenience of users.

# Display help information ###############
display_usage(){
    echo "
    Usage:
           chmod a+x linkFiles.sh  # before the first run
           ./linkFiles.sh [input TSV file]
           ./linkFiles.sh [input TSV file] 1  # Add the second argument '1' to renew existing links.
    The TSV file should not contain a header line. The first column consists of original file paths, and
    the second column consists of link paths:
	    [old name & path]\t[new name & path]\n
    An example of the TSV file:
    ~/data/genome1_1.fasta\t/scratch/input/genome1_unimelb.fna
    ~/data/genome1_2.fasta\t/scratch/input/genome1_zju.fna
    
    Notice a user must ensure the directory is accessible for storing links.
    "
}

if [ -z $1 ]
then
    display_usage
    exit
fi

# Set the override mode ###############
if [ -z $2 ]
then
    override=false
elif [ "$2" -eq "1" ]
then
    override=true
else
    override=false
fi

# Otherwise, make symbolic links following the input file ###############
while read line; do
    if [ ! -z "$line" ]  # Sometimes empty lines are present in the input TSV file, causing an error of ln if keep them untreated.
    then
        IFS=$'\t' read -ra paths <<< "$line"  # split the delimited string into an arrary of two elements
        target="${paths[1]}"
        origin="${paths[0]}"
        if [ ! -L "$target" ]
        then
            ln -s $origin $target
        elif [ "$override" = true ]
        then
            echo "Warning: redirecting existing link $target -> $(readlink ${target}) to ${origin}."
            unlink $target
            ln -s $origin $target
        else
            echo "Warning: skipped existing link $target"
        fi
    fi
done < "$1"  # expect a file name as an input