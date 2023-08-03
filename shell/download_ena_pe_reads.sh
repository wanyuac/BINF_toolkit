#!/bin/bash
# A wrapper for downloading paired-end read sets from the ENA database using ena-file-downloader (github.com/enasequence/ena-ftp-downloader)
# Copyright (C) 2021-2023 Yu Wan <wanyu@microbialsystems.cn>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 29 Mar 2021; latest update: 3 Aug 2023

# Guidance #########################
display_usage() {
    echo "
    Download read files from the ENA database using ena-file-downloader. Please ensure java is accessible in your working environment.
    Command:
        download_ena_pe_reads.sh -d=[path to ena-file-downloader.jar] -t=[a TSV file with two columns: isolate names, read accessions] -o=[Output directory]
    Example commands:
        download_ena_pe_reads.sh -d=\"\$HOME/bin/ena-file-downloader.jar\" -t=\"readsets.tsv\" -o=\"\$PWD\" 1> download_ENA_reads.log 2> download_ENA_reads.err
    "
}

if [ -z "$1" ] || [ "$1" == "-h" ]; then
    display_usage
    exit 0
fi

# Functions ####################
check_dir() {
    if [ ! -d "$1" ]; then
        echo "Create directory $1"
        mkdir -p "$1"
    fi
}

download_reads() {
    p="$1"  # The downloader program
    i="$2"  # Isolate name
    a="$3"  # ENA accession
    tmp_dir="reads_fastq/$a"  # Output directory
    if [ -d "$tmp_dir" ]; then
        echo "Warning: existing temporary directory $tmp_dir is deleted." >&2
        rm -rf "$tmp_dir"
    fi
    java -jar "$p" --accessions="$a" --format=READS_FASTQ --location=$PWD --protocol=FTP --asperaLocation=null  # A directory 'reads_fastq' and a subdirectory "reads_fastq/$j" are created by this command.
    r1="$tmp_dir/${a}_1.fastq.gz"
    r2="$tmp_dir/${a}_2.fastq.gz"
    if [ -f "$r1" ] && [ -f "$r2" ]; then
        mv "$r1" "${i}_1.fastq.gz"
        mv "$r2" "${i}_2.fastq.gz"
        echo "Successfully downloaded paired-end readset of isolate $i (accession: $a)."
    else
        echo "Error: paired-end readset of isolate $i (accession: $a) could not be downloaded." >&2
        echo "Files in directory ${tmp_dir}:" >&2
        ls -1 "${tmp_dir}" >&2
    fi
    rmdir "$tmp_dir"
    sleep 1
}

# Main #########################
# Read arguments
for i in "$@"; do
    case "$i" in
    -t=*)
    accessions="${i#*=}"
    ;;
    -o=*)
    outdir="${i#*=}"
    ;;
    -d=*)
    downloader="${i#*=}"
    ;;
    *)
    ;;
    esac
done

# Check whether the downloader is accessible
if [ ! -f $downloader ]; then
    echo "Error: ena-file-downloader.jar is not accessible at location $downloader" >&2
    exit 1
fi

# Check accession file
if [ ! -f "$accessions" ]; then
    echo "Error: the TSV file of accession numbers is not found." >&2
    exit 1
fi

# Set up the output directory
if [ ! -z "$outdir" ]; then
    check_dir "$outdir"
else
    echo "Error: $outdir was not found." >&2
    exit 1
fi

# Download reads
cd "$outdir"

while read line; do  # Read through the input TSV file line-by-line
    if [ ! -z "$line" ]; then
        IFS=$'\t' read -r -a fields <<< "$line"
        download_reads "$downloader" "${fields[0]}" "${fields[1]}"  # Downloader, isolate name, ENA accession
    fi
    rmdir reads_fastq  # A directory created by the downloader; directory 'logs' (also created by the downloader) is left in the output directory.
done < "$accessions"  # Expect a file name as an input
