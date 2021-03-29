#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyu@microbialsystems.cn>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 29/3/2021

# Guidance #########################
display_usage() {
    echo "
    Download read files from the ENA via FTP.
    Input: a three-column TSV file: [isolate name]\t[FTP address of read set 1]\t[FTP address of read set 2]
    Example commands:
        cd [read directory]
        download_ena_reads.sh readsets.tsv 1> download_ENA_reads.log 2> download_ENA_reads.err  # Reads will be downloaded to \$PWD.
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Main #########################
while read line; do
    IFS=$'\t' read -r -a fields <<< "$line"
    readfile1="${fields[0]}_1.fastq.gz"
    readfile2="${fields[0]}_2.fastq.gz"
    echo -e "Download ${fields[1]} and save it as ${readfile1}."
    wget -O ${readfile1} ${fields[1]}
    sleep 1
    echo -e "Download ${fields[2]} and save it as ${readfile2}."
    wget -O ${readfile2} ${fields[2]}
    sleep 1
done < "$1"  # expect a file name as an input
