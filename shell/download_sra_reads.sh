#!/bin/bash
# Copyright (C) 2020-2025 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 11 March 2020; last modification: 17 April 2025
# Important update on 16/4/2025: added fastq-dump arguments "--skip-technical --clip --dumpbase --read-filter pass"
# according to https://edwards.flinders.edu.au/fastq-dump/. (Thanks to Sophie Mannix for pointing this out)

# Help information #########################
show_help() {
    echo "
    Download files of paired-end short reads from the NCBI SRA database.
    Arguments:
        -d: (optional) Directory that contains the program fastq-dump.
        -o: Output directory (no forward slash).
        -a: A comma-delimited string of target accession numbers (SRR*)
        -f: A single-column text file of SRR numbers or a two-column CSV file of genome names (1st column)
            and SRR numbers (2nd column).
        -r: A logical flag turning on replacement of SRR numbers with genome names for read files.
        -s: A logical flag notifying this script that the reads to be downloaded are single-end.
        -u: Skip the dos2unix step if your input file is known to follow the Unix-style line ending.
    Example command:
        ./download_sra_reads.sh -d=\"\$HOME/bin/sra_toolkit/bin\" -o=\"\$PWD\" -f=readsets.csv -r
    Note that:
        1. The -a argument is ignored when the -f argument is set.
        2. Newline characters in the input file must be '\n' rather than '\r\n'.
        3. Fastq-dump sometimes fails in downloading or parsing read files. Remember to check the STDOUT output
           of download_sra.sh after each run.
        4. Dependency: SRA Toolkit v3.0.6 and later versions, dos2unix
    "
}

if [[ $# -eq 0 ]]; then  # When $1 does not exist (Error "$1: unbound variable" arises when use `if [[ $# -eq 0 ]] || [ "$1" = "-h" ]`). Do not use `-z "$1"`.
    show_help
    exit 0
elif [ "$1" = "-h" ]; then
    show_help
    exit 0
fi

# Remove leading and trailing whitespace (spaces, tabs, etc.) while handling carriage returns (\r) and newlines (\n) #########################
function trim_whitespace() {
    echo -e "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}

# Main function #########################
out_dir='.'  # Output directory
replace_names=false  # By default, do not replace accession numbers with genome names.
paired_end=true  # Assumes all read files are paired-end.
read_file=false  # Assumes that by default accessions are not provided in a file.
unix_format=false  # Assumes the input file has non-Unix-style line endings ('\n\r' etc).

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
        -d=*)
        program_dir="${i#*=}"
        ;;
        -o=*)
        out_dir="${i#*=}"
        ;;
        -u)
        unix_format=true
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

# Check the output directory
if [ ! -d "$out_dir" ]; then
    echo "Making the output directory: $out_dir"
    mkdir -p "$out_dir"
fi

# Load module
if [ ! -z "$program_dir" ]; then
    export PATH="${program_dir}:$PATH"
fi

# Download and parse read files
error_count=0
successes=0
if [ "$read_file" = true ]; then
    if [ "$unix_format" = false ]; then
        dos2unix "$acc_list"
    fi
        
    # Read lines of the input file into an array, skipping empty lines
    # https://stackoverflow.com/questions/15685736/how-to-extract-a-particular-element-from-an-array-in-bash
    echo "Reading accession numbers from file $acc_list"
    mapfile -t lines_array < "$acc_list"
    echo "Read ${#lines_array[@]} entries from file $acc_list"
    
    if [ "$replace_names" = true ]; then
        echo -e "Accession numbers in file names will be replaced by genome names.\n"
        for line in "${lines_array[@]}"; do
            echo "Parsing line '${line}'"
            IFS=',' read -r -a line_fields <<< "${line}"
            genome="${line_fields[0]}"
            accession="${line_fields[1]}"
            genome="$(trim_whitespace "$genome")"
            accession="$(trim_whitespace "$accession")"
            echo "Accession: ${accession}; Genome: ${genome}"
            if [ "$paired_end" = true ]; then  # Paired-end reads
                echo "Download $accession and rename files as ${genome}_1.fastg.gz and ${genome}_2.fastq.gz."
                fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession"  # Download and split the read file, and create the output directory if necessary
                f1="${out_dir}/${accession}_pass_1.fastq"
                f2="${out_dir}/${accession}_pass_2.fastq"
                if [ -f "$f1" ] && [ -f "$f2" ]; then
                    echo "Compress $f1 and $f2 and rename output files"
                    gzip "$f1"
                    gzip "$f2"
                    mv "${f1}.gz" "${out_dir}/${genome}_1.fastq.gz"
                    mv "${f2}.gz" "${out_dir}/${genome}_2.fastq.gz"
                    echo "Completed the compression and name changes for $accession"
                    ((successes++))
                else
                    echo "Error: $f1 and/or $f2 could not be downloaded." >&2
                    ((error_count++))
                fi
            else  # Single-end reads
                echo "Download ${accession} and rename it as ${genome}.fastq.gz."
                fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession"
                f1="${out_dir}/${accession}_pass.fastq"
                if [ -f "$f1" ]; then
                    gzip "$f1"
                    mv "${f1}.gz" "${out_dir}/${genome}.fastq.gz"
                    ((successes++))
                else
                    echo "Error: $f1 could not be downloaded." >&2
                    ((error_count++))
                fi
            fi
            echo -e "Finished processing ${accession}.\n"
            sleep 2  # Pause, to avoid too many connection requests to NCBI's server.
        done
    else  # A single-column input file
        echo -e "Use accession numbers as filenames of downloaded read sets.\n"
        for accession in "${lines_array[@]}"; do
            accession="$(trim_whitespace "$accession")"
            echo "Downloading read set ${accession}."
            fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir ${out_dir} --split-3 ${accession}
            if [ "${paired_end}" = true ]; then
                f1="${out_dir}/${accession}_pass_1.fastq"
                f2="${out_dir}/${accession}_pass_2.fastq"
                if [ -f "$f1" ] && [ -f "$f2" ]; then
                    gzip "$f1"
                    gzip "$f2"
                    ((successes++))
                else
                    echo "Error: $f1 and/or $f2 could not be downloaded." >&2
                    ((error_count++))
                fi
            else
                f1="${out_dir}/${accession}_pass.fastq"
                if [ -f "$f1" ]; then
                    gzip "$f1"
                    ((successes++))
                else
                    echo "Error: $f1 could not be downloaded." >&2
                    ((error_count++))
                fi
            fi
            echo -e "Finished processing ${accession}.\n"
            sleep 2
        done
    fi
else  # When accession numbers come from the -a parameter
    echo -e "Reading accession numbers from the parameter '-a'.\n"
    for accession in "${accessions[@]}"; do  # Filename replacement is not supported under this mode.
        accession="$(trim_whitespace "$accession")"
        echo "Downloading and parsing ${accession}."
        fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir ${out_dir} --split-3 $accession
        if [ "${paired_end}" = true ]; then
            f1="${out_dir}/${accession}_1.fastq"
            f2="${out_dir}/${accession}_2.fastq"
            if [ -f "$f1" ] && [ -f "$f2" ]; then
                gzip "$f1"
                gzip "$f2"
            else
                echo "Error: $f1 and/or $f2 could not be downloaded." >&2
            fi
        else
            f1="${out_dir}/${accession}.fastq"
            if [ -f "$f1" ]; then
                gzip "$f1"
            else
                echo "Error: $f1 could not be downloaded." >&2
            fi
        fi
        echo -e "Finished processing ${accession}.\n"
        sleep 2
    done
fi

echo "Finished all download tasks. Downloaded $successes readsets and failed to download $error_count readsets"
exit 0
