#!/bin/bash
# This script concatenates assembly statistics from our lab's previous assembly pipeline.
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# Creation: 31/10/2016; the latest version: 28/11/2017

display_usage(){
    echo "
    Concatenates *_contigStats.txt files that are generated using our lab's contigMetrics.py.
    Usage: bash catContigStats.sh ./assemblies/
    Outputs:
        contigStats_files.txt
        contigStats_combined.csv
    "
}

HEADER="contigFile,numContigs,totalBases,N50,smallest,lowerQ,median,upperQ,largest"
FILE_LIST="contigStats_files.txt"
STATS="contigStats_combined.csv"

# Check argument ##########
if [ -z $1 ]; then
    echo "Error: a subject directory must be provided."
    display_usage
    exit
fi

# Find all contigStats.txt files ##########
find $1 -name *_contigStats.txt -type f > $FILE_LIST
echo "There are `cat ${FILE_LIST} | wc -l` genomes."

# Print contig statistics into a CSV file ##########
echo $HEADER > $STATS

# Extract the second line of every file and appends it to the CSV file ==========
files=`cat ${FILE_LIST}`
n=0

for f in ${files}; do
    r=`cat ${f} | wc -l`
    if [ "$r" -eq "2" ]; then
        tail -n 1 $f >> $STATS
        ((n++))
    else
        echo "Warning: ${f} does not contain contig statistics."
    fi
done

echo "Success: ${n} lines of statistics have been transferred into ${STATS}."
