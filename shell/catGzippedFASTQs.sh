#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First edition: 25 Aug 2021; the latest update: 25 Aug 2021

# User guide ####################
display_usage() {
    echo "
    Concatenate *.fastq.gz in each subdirectory into a single fastq.gz file. Useful for concatenating
    Guppy's output sequence files.
    Usage:
        bash catGzippedFASTQs.sh [input parental directory] [output directory] [inputs.tsv]
    For example: ./catGzippedFASTQs.sh barcodes.tsv ~/fastq/pass ~/fastq/concat > cat_fastqs.log
    Input TSV file inputs.tsv consists of two columns: subdirectory name, isolate name (for the output file).
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Main utility ####################
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

n=0  # Count the number of subdirectories visited
while read line  # Please ensure every line in the input TSV file is ended with a newline character.
do
    if [ ! -z "$line" ]
    then
        IFS=$'\t' read -r -a fields <<< "$line"  # Parse the line into two fields by '\t'.
        input_subdir="$indir/${fields[0]}"
        if [ ! -d "$input_subdir" ]
        then
            echo "Skip inaccessible input directory $input_subdir" >&2
        else
            output="$outdir/${fields[1]}.fastq.gz"
            k=$(ls -1 $input_subdir/*.fastq.gz | wc -l)
            echo "Concatenate $k .fastq.gz files from $input_subdir into $output"
            zcat $input_subdir/*.fastq.gz | gzip > $output  # Slower than 'cat *.fastq.gz' but generates a smaller file.
            (( n++ ))
        fi
    fi
done < "$3"

echo "FASTQ files of $n samples have been successfully concatenated."