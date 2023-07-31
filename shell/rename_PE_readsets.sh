#!/bin/bash
# Renaming Illumina paired-end readsets via symbolic links.
# Copyright (C) 2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 31 July 2023; latest update: 31 July 2023

display_usage() {
    echo "Rename paired-end readsets via symbolic links.
    Command line:
    rename_PE_readsets.sh [mapping file] [directory of original readsets] [directory for links]
    The mapping file is TSV-delimited and does not have any header. It consists of two columns:
    original name and new name, respectively. This script assumes filenames of readsets have
    suffices _1.fastq.gz and _2.fastq.gz"
}

if [ -z "$1" ] || [ $1 = "-h" ]; then
    display_usage
    exit
fi

dir_in="$2"
dir_out="$3"

if [ ! -d "$dir_in" ]; then
    echo "Error: input directory $dir_in was not found." >&2
    exit
fi

if [ ! -d "$dir_out" ]; then
    echo "Create output directory $dir_out"
    mkdir -p "$dir_out"
fi

while read -r line; do
    IFS=$'\t' read -r -a vals <<< "$line"
    i="${vals[0]}"  # Original name
    j="${vals[1]}"  # New name
    r1="$dir_in/${i}_1.fastq.gz"
    r2="$dir_in/${i}_2.fastq.gz"
    if [ -f "$r1" ] && [ -f "$r2" ]; then
        t1="$dir_out/${j}_1.fastq.gz"
        t2="$dir_out/${j}_2.fastq.gz"
        echo -e "$r1 -> $t1\t$r2 -> $t2"
        ln -s "$r1" "$t1"
        ln -s "$r2" "$t2"
    else
        echo "Error: $r1 or $r2 were not accessible. No links were created for sample ${i}." >&2
    fi
done < "$1"
