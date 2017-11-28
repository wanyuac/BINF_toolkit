#!/bin/bash

# Concatenates files of assembly statistics that are generated using our lab's script contigMetrics.py.
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# Development history: 31/10/2016, 20/11/2017
# Previous name: catContigStats.sh

display_usage(){
    echo "
    Concatenates files of assembly statistics that are generated using our lab's script contigMetrics.py.
    Usage: bash catAssemblyStats.sh ./assemblies/
    Outputs:
        assemblyStats_files.txt
        assemblyStats.tsv
    "
}

#Constants for the previous version of contigMetrics.py
#HEADER="contigFile,numContigs,totalBases,N50,smallest,lowerQ,median,upperQ,largest"  # previous version
#STATS="contigStats_combined.csv"
#FILE_LIST="contigStats_files.txt"
#
#Previous outputs:
#    contigStats_files.txt
#    contigStats_combined.csv

HEADER="Assembly\tContig_number\tN50\tQ1\tQ2\tQ3\tMean\tSmallest\tLargest\tLength"
STATS="assemblyStats.tsv"
FILE_LIST="assemblyStats_files.txt"
NAME_PATTERN="assembly_stats.tsv"  # previous name: *_contigStats.txt

# Check argument ##########
if [ -z $1 ]; then
    echo "Error: a subject directory must be provided."
    display_usage
    exit
fi

# Find all contigStats.txt files ##########
find $1 -name $NAME_PATTERN -type f > $FILE_LIST
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
