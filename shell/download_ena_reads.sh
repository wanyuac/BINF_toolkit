#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyu@microbialsystems.cn>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 29/3/2021; latest update: 31/3/2021

# Guidance #########################
display_usage() {
    echo "
    Download read files from the ENA via FTP (method: wget), skipping downloaded files in the local folder.
    Input: a three-column TSV file: [isolate name]\t[FTP address of read set 1]\t[FTP address of read set 2]
    Command format:
        download_ena_reads.sh [TSV file for URLs] [Output directory]
    Example commands:
        download_ena_reads.sh readsets.tsv ./reads 1> download_ENA_reads.log 2> download_ENA_reads.err  # Reads will be downloaded to \$PWD.
    "
}

# Construct the downloader function download_and_check ####################
display_cmd() {
    # A subordinate function of download_reads.
    echo "
# Download $2 ==> $1:
wget --timeout=20 --waitretry=2 --tries=3 --retry-connrefused --no-verbose --no-glob --output-document $1 $2
sleep 6
    "
}

download_reads() {
    # $1: output filename; $2: FTP address of the remote read file
    # A subordinate function of download_and_check.
    if [ ! -f "$1" ]; then
        display_cmd $1 $2
        wget --timeout=20 --waitretry=2 --tries=3 --retry-connrefused --no-verbose --no-glob --output-document $1 $2
        sleep 6  # Waiting for 2-4 seconds may not be able to avoid the problem of "Gateway timeout"
    else
        echo "# Skip downloading for $1 as it exists."
    fi
}

download_and_check() {
    # $1: output filename; $2: FTP address of the remote read file
    # retries=`seq 1 1 3`
    download_reads $1 $2  # Rety three times if the problem of "Gateway timeout" keeps occurring.
    #for i in ${retries[@]}; then
    for i in 1 2 3; do
        if [ -s "$1" ]; then
            break
        else
            rm -f $1
            download_reads $1 $2  # Download the file again if the existing one has a size of zero (due to the error of "Gateway timeout")
        fi
    done
}

# Main #########################
if [ -z $1 ]; then
    display_usage
    exit
fi

# Set up output directory
if [ -z $2 ]; then
    outdir=$PWD
else
    outdir=$2
fi

if [ ! -d "$outdir" ]; then
    echo "# Create read directory ${outdir}."
    mkdir $outdir
else
    echo "# Read directory $outdir exists."
fi

while read line; do  # Read through the input TSV file line-by-line
    IFS=$'\t' read -r -a fields <<< "$line"
    isolate=${fields[0]}
    url1=${fields[1]}
    url2=${fields[2]}
    readfile1="${outdir}/${isolate}_1.fastq.gz"
    readfile2="${outdir}/${isolate}_2.fastq.gz"
    download_and_check $readfile1 $url1
    download_and_check $readfile2 $url2
done < "$1"  # Expect a file name as an input
