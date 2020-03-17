#!/bin/bash
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 11 March 2020
# Last modification: 17 March 2020

# Help information #########################
show_help() {
    echo "
    Download files of paired-end short reads from the NCBI SRA database.
    Arguments:
        -m: (optional) Name of an sra-toolkit module, which contains the program fastq-dump.
        -o: Output directory (no forward slash).
        -a: A comma-delimited string of target accession numbers (SRR*)
        -f: A single-column text file of SRR numbers or a two-column tab-delimited file of SRR numbers (1st column)
            and genome names (2nd column).
        -r: A logical flag turning on replacement of SRR numbers with genome names for read files.
        -s: A logical flag notifying this script that the reads to be downloaded are single-end.
    Example command:
        download_sra.sh -m='sra-toolkit/2.8.1-3' -o='reads' -a='SRR00001,SRR00002,SRR00003' > download.log
        download_sra.sh -m='sra-toolkit/2.8.1-3' -o='reads' -f='accessions.tsv' -r > download.log
    Note that:
        1. The -a argument is ignored when the -f argument is set.
        2. Newline characters in the input file must be '\n' rather than '\r\n'.
        3. Fastq-dump sometimes fails in downloading or parsing read files. Remember to check the STDOUT output
           of download_sra.sh after each run.
    "
}

if [ -z "$1" ] || [ $1 = "-h" ]; then
    show_help
    exit
fi

# Main function #########################
env_module=''  # Name of sra-toolkit
out_dir='.'  # Output directory
replace_names=false  # By default, do not replace accession numbers with genome names.
paired_end=true  # Assumes all read files are paired-end.
read_file=false  # Assumes that by default accessions are not provided in a file.

# Read arguments
for i in "$@"; do
    case $i in
        -a=*)
        accessions="${i#*=}"
        accessions=( `echo $accessions | tr "," "\n"` )
        ;;
        -f=*)
        acc_list="${i#*=}"  # A file listing accession numbers and probably genome names as well
        read_file=true
        ;;
        -r)
        replace_names=true
        ;;
        -s)
        paired_end=false  # Single-end reads
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

# Check the output directory
if [ ! -d "${out_dir}" ]; then
    echo "Making the output directory: ${out_dir}"
    mkdir -p ${out_dir}
fi

# Load module
if [ ! -z "${env_module}" ]; then
    echo "Loading environment module ${env_module}."
    module load ${env_module}
fi

# Download and parse read files
if [ "${read_file}" = true ]; then
    echo "Reading accession numbers from file ${acc_list}."
    
    # Read lines of the input file into an array, skipping empty lines
    # https://stackoverflow.com/questions/15685736/how-to-extract-a-particular-element-from-an-array-in-bash
    IFS=$'\n' read -d '' -r -a lines_array < "${acc_list}"
    
    if [ "${replace_names}" = true ]; then
        echo -e "Genome names will substitute for accession numbers for read sets.\n"
        for line in "${lines_array[@]}"; do
            IFS=$'\t' read -r -a line_fields <<< "${line}"
            accession="${line_fields[0]}"
            genome="${line_fields[1]}"
            if [ "${paired_end}" = true ]; then  # Paired-end reads
                echo "Downloading ${accession} and rename files as ${genome}_1.fastg.gz and ${genome}_2.fastq.gz."
                fastq-dump --readids --outdir ${out_dir} --gzip --split-3 ${accession}  # Download and split the read file
                mv ${out_dir}/${accession}_1.fastq.gz ${out_dir}/${genome}_1.fastq.gz
                mv ${out_dir}/${accession}_2.fastq.gz ${out_dir}/${genome}_2.fastq.gz
            else  # Single-end reads
                echo "Download ${accession} and rename it as ${genome}.fastq.gz."
                fastq-dump --readids --outdir ${out_dir} --gzip --split-3 ${accession}
                mv ${out_dir}/${accession}.fastq.gz ${out_dir}/${genome}.fastq.gz
            fi
            echo -e "Finished processing ${accession}.\n"
            sleep 2  # Pause, to avoid too many connection requests to NCBI's server.
        done
    else  # A single-column input file
        echo -e "Use accession numbers as filenames of downloaded read sets.\n"
        for accession in "${lines_array[@]}"; do
            echo "Downloading read set ${accession}."
            fastq-dump --readids --outdir ${out_dir} --gzip --split-3 ${accession}
            echo -e "Finished processing ${accession}.\n"
            sleep 2
        done
    fi
else  # When accession numbers come from the -a parameter
    echo -e "Reading accession numbers from the parameter '-a'.\n"
    for i in "${accessions[@]}"; do  # Filename replacement is not supported under this mode.
        echo "Downloading and parsing ${i}."
        fastq-dump --readids --outdir ${out_dir} --gzip --split-3 $i
        echo -e "Finished processing ${i}.\n"
        sleep 2
    done
fi

echo 'All download tasks have been finished successfully.'
