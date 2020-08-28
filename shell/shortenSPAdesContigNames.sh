#!/bin/bash
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Release: 28 Aug 2020

display_usage(){
    echo "
    Shorten contig names from SPAdes by substituting white spaces for '_length_' in sequence headers,
    so the latter part of sequence description will be ignored by Prokka and some other software.
    This is particularly useful for Prokka annotation because a long contig name consumes all space
    between the contig name and length in the output GenBank file, causing a problem to SniEff and so on.
    
    Command line: bash shortenSPAdesContigNames.sh [input FASTA file] > [new FASTA file]
    "
}

if [ -z $1 ]; then
    display_usage
    exit
else
    sed 's/_length_/ /g' $1
fi
