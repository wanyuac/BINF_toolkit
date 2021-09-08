#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First edition: 8 Sep 2021; the latest update: 8 Sep 2021
# This script is derived from catGzippedFASTQsPerDirectory.sh.

# User guide ####################
display_usage() {
    echo "
    Concatenate *.fastq.gz of each sample into a single fastq.gz file.
    Usage:
        bash catGzippedFASTQsPerSample.sh [input parental directory] [output directory] [a list of input sample names]
    For example: ./catGzippedFASTQsPerSample.sh ~/fastq/pass isolates.txt ~/fastq/concat &> cat_fastqs.log
    There is one sample name per line in the input sample-name list.
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Main utility ####################

# 1. Set up directories ===============
indir="$1"
outdir="$2"

if [ ! -d "$indir" ]
then
    echo "Error: input parental directory $indir does not exist." >&2  # Print to standard error
    exit
fi

if [ ! -d "$outdir" ]
then
    echo "Making output directory $outdir"
    mkdir $outdir
fi

# 2. Concatenate read files ===============
n=0  # The counter of samples processed
while read i  # Please ensure every line in the input TSV file is ended with a newline character.
do
    if [ ! -z "$i" ]  # Skip empty lines
    then
        # Users may customise the following two commands to match their filenames.
        ra="$indir/*_${i}A-[1,2].bacterial-fastq-only.ngsservice.processed.R"  # There should be only a single match.
        rb="$indir/*_${i}B-[1,2].bacterial-fastq-only.ngsservice.processed.R"  # The same as above.
        ra1=`ls -1 ${ra}1.fastq.gz`
        ra2=`ls -1 ${ra}2.fastq.gz`
        rb1=`ls -1 ${rb}1.fastq.gz`
        rb2=`ls -1 ${rb}2.fastq.gz`

        # Concatenate the read files of the current sample
        if [ -f "$ra1" ] && [ -f "$ra2" ] && [ -f "$rb1" ] && [ -f "$rb2" ]
        then
            echo "Process read files of isolate $i"
            echo "    Concatenating $ra1 and $rb1"
            zcat $ra1 $rb1 | gzip > $outdir/${i}_1.fastq.gz  # Slower than 'cat *.fastq.gz' but generates a smaller file.
            echo "    Concatenating $ra2 and $rb2"
            zcat $ra2 $rb2 | gzip > $outdir/${i}_2.fastq.gz
            (( n++ ))
        else
            echo "Skip file concatenation for sample $i due to absence of one or more read files." >&2
        fi
    fi
done < "$3"  # Read sample names one-by-one from the input list

echo "FASTQ files of $n samples have been successfully concatenated."