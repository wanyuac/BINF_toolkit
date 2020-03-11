#!/bin/bash
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 11 March 2020

# Help information #########################
show_help() {
    echo "
    Download files of paired-end short reads from the NCBI SRA database.
    Arguments:
        -m: (optional) Name of an sra-toolkit module, which contains the program fastq-dump.
        -o: Output directory (no forward slash).
        -a: A comma-delimited string of target accession numbers (SRR*)
    Example command:
        download_sra.sh -m='sra-toolkit/2.8.1-3' -o='reads' -a='SRR00001,SRR00002,SRR00003'
    "
}

if [ -z "$1" ] || [[ $1 == -h ]]; then
	show_help
	exit
fi

# Main function #########################
env_module=''  # Name of sra-toolkit
out_dir='.'  # Output directory

# Read arguments
for i in "$@"; do
    case $i in
        -a=*)
        accessions="${i#*=}"
        accessions=( `echo $accessions | tr "," "\n"` )
        ;;
		-m=*)
        env_module="${i#*=}"
		;;
        -o=*)
        out_dir="${i#*=}"
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

# Load module
if [ ! -z "${env_module}" ]; then
	echo "Loading environment module ${env_module}."
	module load ${env_module}
fi

# Download and parse read files
cd ${out_dir}

for i in "${accessions[@]}"; do
	echo "Downloading ${i}."
	fastq-dump --gzip --readids --split-3 $i
done

echo 'All download tasks have been finished successfully.'
