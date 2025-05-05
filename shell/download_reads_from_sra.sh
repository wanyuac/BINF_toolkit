#!/bin/bash
# Copyright (C) 2020-2025 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 11 March 2020; last modification: 5 May 2025
# Important update on 16/4/2025: added fastq-dump arguments "--skip-technical --clip --dumpbase --read-filter pass"
# according to https://edwards.flinders.edu.au/fastq-dump/. (Thanks to Sophie Mannix for pointing this out)

# Help information #########################
show_help() {
    echo "
    Download files of sequencing reads from the NCBI SRA database.
    Arguments:
        -d: (optional) Directory that contains the program fastq-dump (does not need to end by a slash character '/').
        -o: Output directory (no forward slash). Default: ${HOME}/SRA_reads.
        -a: A comma-delimited string of target accession numbers (SRR*)
        -f: A single-column text file of SRR numbers or a two-column CSV file of genome names (1st column)
            and SRR numbers (2nd column).
        -r: A logical flag turning on replacement of SRR numbers with genome names for read files.
        -s: A logical flag notifying this script that the reads to be downloaded are single-end.
        -u: Skip the dos2unix step if your input file is known to follow the Unix-style line ending.
        -p: Prefix of the log file (Markdown format) in the output directory (default: download_reads_from_sra_[date (YYYY-MM-DD)]_[HH-MM-SS])
    Example command:
        ./download_reads_from_sra.sh -d=\"\$HOME/bin/sra_toolkit/bin\" -o=\"\$PWD\" -f=readsets.csv -r -p=download_reads_from_sra
    Note that:
        1. The -a argument is ignored when the -f argument is set.
        2. Newline characters in the input file must be '\n' rather than '\r\n'.
        3. Fastq-dump sometimes fails in downloading or parsing read files. Remember to check the log and error files after each run.
        4. Please ensure genome names are unique throughout your dataset. Otherwise, files of the same names may be overridden at the renaming step.
        5. Dependency: SRA Toolkit v3.0.6 and later versions, dos2unix.
        6. This script was called download_sra_reads.sh.
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

function time_stamp() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]"
}

function write_log() {
    echo "$1" >> "$log"  # The global variable log will be defined later.
}

# Main function #########################
out_dir="${HOME}/SRA_reads"  # The default output directory
replace_names=false  # By default, do not replace accession numbers with genome names.
paired_end=true  # Assumes all read files are paired-end.
read_file=false  # Assumes that by default accessions are not provided in a file.
unix_format=false  # Assumes the input file has non-Unix-style line endings ('\n\r' etc).
prefix=download_reads_from_sra
suffix=$(date +"%Y-%m-%d_%H-%M-%S")  # Name suffix of log and error files
wait_time=2  # In seconds. Time to pause before the start of a new download iteration

# Read arguments
for i in "$@"; do
    case $i in
        -a=*)
        accessions="${i#*=}"
        accessions=( $(echo $accessions | tr ',' '\n') )
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
        -p=*)
        prefix="${i#*=}"
        ;;
        -w=*)
        wait_time="${i#*=}"
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

log="${out_dir}/${prefix}_${suffix}.log"

# Check the output directory
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
    write_log "# Global configurations"
    write_log "$(time_stamp) Created the output directory: ${out_dir}."
fi

# Load module
if [ ! -z "$program_dir" ]; then
    export PATH="${program_dir}:$PATH"
fi

# Download and parse read files
write_log "Waiting time between consecutive download iterations: $wait_time seconds."
accession_count=0  # Assign a default value to accession_count in case neither -f nor -a is set.
successes=0
failures=0

if [ "$read_file" = true ]; then
    if [ ! -f "$acc_list" ]; then
        write_log "$(time_stamp) Error: $acc_list was not found."
        exit 1
    fi
    if [ "$unix_format" = false ]; then
        dos2unix "$acc_list" >> "$log" 2>&1
    fi
        
    # Read lines of the input file into an array, skipping empty lines
    # https://stackoverflow.com/questions/15685736/how-to-extract-a-particular-element-from-an-array-in-bash
    mapfile -t lines_array < "$acc_list"
    accession_count="${#lines_array[@]}"  # Number of accessions
    write_log "$(time_stamp) Imported $accession_count entries from file ${acc_list}."
    if [ "$replace_names" = true ]; then
        write_log "Run accessions in FASTQ filenames will be replaced by genome names."
        echo >> "$log"  # Create an empty line in the log file
        write_log "# Task records"
        for line in "${lines_array[@]}"; do
            # write_log "$(time_stamp) Parse line '${line}'"  # This command is used for debugging
            IFS=',' read -r -a line_fields <<< "${line}"
            genome="${line_fields[0]}"  # Genome or isolate name
            accession="${line_fields[1]}"  # Run accession
            genome="$(trim_whitespace "$genome")"  # Sometimes people include whitespaces or tab characters in some genome names or accessions by accident.
            accession="$(trim_whitespace "$accession")"
            if [[ -z "$genome" || -z "$accession" ]]; then
                write_log "$(time_stamp) Warning: Skipping malformed line: '$line'"
                continue
            fi
            fastq_original_prefix="${out_dir}/${accession}"
            fastq_renamed_prefix="${out_dir}/${genome}"
            # write_log "$(time_stamp)  genome name ${genome} and its Run accession ${accession}"  # For debugging
            if [ "$paired_end" = true ]; then  # Paired-end mode
                write_log "## Download paired-end reads under Run accession $accession of genome ${genome}"
                write_log "$(time_stamp) Run fastq-dump to download the reads."
                if fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession" >> "$log" 2>&1; then  # Download and split the read file, and create the output directory if necessary
                    write_log "$(time_stamp) fastq-dump command finished for $accession of genome ${genome}."
                else
                    write_log "$(time_stamp) fastq-dump failed for $accession."
                    ((failures++))
                    continue
                fi
                f1="${fastq_original_prefix}_pass_1.fastq"
                f2="${fastq_original_prefix}_pass_2.fastq"
                f1_renamed="${fastq_renamed_prefix}_1.fastq"
                f2_renamed="${fastq_renamed_prefix}_2.fastq"
                if [ -f "$f1" ] && [ -f "$f2" ]; then
                    #write_log "$(time_stamp) Compress FASTQ files $f1 and $f2"
                    #gzip "$f1"  # Removed the gzip step because occasionally it fails to execute on my university's HPC
                    #gzip "$f2"
                    #if [ -f "${f1}.gz" ] && [ -f "${f2}.gz"]; then
                        #mv "${f1}.gz" "$f1_renamed"
                        #mv "${f2}.gz" "$f2_renamed"
                        #write_log "$(time_stamp) Successfully compressed and renamed FASTQ files: ${f1}.gz -> ${f1_renamed}; ${f2}.gz -> ${f2_renamed}"
                    mv "$f1" "$f1_renamed"
                    mv "$f2" "$f2_renamed"
                    if [ -f "$f1_renamed" ] && [ -f "$f2_renamed" ]; then
                        write_log "$(time_stamp) Successfully renamed FASTQ files: $f1 -> ${f1_renamed}; $f2 -> ${f2_renamed}."
                        ((successes++))
                    else
                        write_log "$(time_stamp) Error: FASTQ files $f1 and/or $f2 could not be renamed."
                        ((failures++))
                    fi
                else
                    write_log "$(time_stamp) Error: $f1 and/or $f2 could not be created by fastq-dump. Skip the step of renaming FASTQ files."
                    ((failures++))
                fi
            else  # Single-end mode (for instance, PacBio or Nanopore reads)
                write_log "## Download single-end read set $accession of genome ${genome}"
                write_log "$(time_stamp) Run fastq-dump to download the reads."
                if fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession" >> "$log" 2>&1; then
                    write_log "$(time_stamp) fastq-dump command finished for $accession of genome ${genome}."
                else
                    write_log "$(time_stamp) fastq-dump failed for $accession."
                    ((failures++))
                    continue
                fi
                f1="${fastq_original_prefix}_pass.fastq"
                f1_renamed="${fastq_renamed_prefix}.fastq"
                if [ -f "$f1" ]; then
                    #gzip "$f1"
                    mv "$f1" "$f1_renamed"
                    if [ -f "$f1_renamed" ]; then
                        write_log "$(time_stamp) Successfully renamed FASTQ file: $f1 -> ${f1_renamed}."
                        ((successes++))
                    else
                        write_log "$(time_stamp) Error: FASTQ file $f1 could not be renamed."
                        ((failures++))
                    fi
                else
                    write_log "$(time_stamp) Error: $f1 could not be created by fastq-dump."
                    ((failures++))
                fi
            fi
            echo >> "$log"  # Add an empty line before the start of a new download iteration
            sleep "$wait_time"  # Pause, to avoid too many connection requests to NCBI's server.
        done
    else  # A single-column input file of Run accessions
        write_log "Use Run accessions as filenames of downloaded read sets."
        echo >> "$log"
        write_log "# Task records"
        for line in "${lines_array[@]}"; do
            accession="$(trim_whitespace "$line")"
            if [ -z "$accession" ]; then
                write_log "$(time_stamp) Warning: Skipping malformed line: '$line'"
                continue
            fi
            fastq_original_prefix="${out_dir}/${accession}"
            write_log "## Download reads under Run accession $accession"
            write_log "$(time_stamp) Run fastq-dump to download reads under ${accession}."
            if fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession" >> "$log" 2>&1; then
                write_log "$(time_stamp) fastq-dump command finished for ${accession}."
            else
                write_log "$(time_stamp) fastq-dump failed for $accession."
                ((failures++))
                continue
            fi
            if [ "$paired_end" = true ]; then  # Paired-end reads
                f1="${fastq_original_prefix}_pass_1.fastq"
                f2="${fastq_original_prefix}_pass_2.fastq"
                if [ -f "$f1" ] && [ -f "$f2" ]; then
                    #gzip "$f1"
                    #gzip "$f2"
                    write_log "$(time_stamp) Successfully created FASTQ files $f1 and ${f2}."
                    ((successes++))
                else
                    write_log "$(time_stamp) Error: $f1 and/or $f2 could not be created by fastq-dump."
                    ((failures++))
                fi
            else  # Single-end reads
                f1="${fastq_original_prefix}_pass.fastq"
                if [ -f "$f1" ]; then
                    #gzip "$f1"
                    write_log "$(time_stamp) Successfully created FASTQ file ${f1}."
                    ((successes++))
                else
                    write_log "$(time_stamp) Error: $f1 could not be created by fastq-dump."
                    ((failures++))
                fi
            fi
            echo >> "$log"
            sleep "$wait_time"
        done
    fi
else  # When accession numbers come from the -a parameter
    accession_count="${#accessions[@]}"
    write_log "$(time_stamp) Imported $accession_count entries from the -a argument."
    echo >> "$log"
    write_log "# Task records"
    for entry in "${accessions[@]}"; do  # Filename replacement is not supported under this mode.
        accession="$(trim_whitespace "$entry")"
        if [ -z "$accession" ]; then
            write_log "$(time_stamp) Warning: Skipping malformed entry: '$entry'"
            continue
        fi
        fastq_original_prefix="${out_dir}/${accession}"
        write_log "## Download reads under Run accession $accession"
        if fastq-dump --readids --skip-technical --clip --dumpbase --read-filter pass --outdir "$out_dir" --split-3 "$accession" >> "$log" 2>&1; then
            write_log "$(time_stamp) fastq-dump command finished for ${accession}."
        else
            write_log "$(time_stamp) fastq-dump failed for $accession."
            ((failures++))
            continue
        fi
        if [ "${paired_end}" = true ]; then
            f1="${fastq_original_prefix}_1.fastq"
            f2="${fastq_original_prefix}_2.fastq"
            if [ -f "$f1" ] && [ -f "$f2" ]; then
                #gzip "$f1"
                #gzip "$f2"
                write_log "$(time_stamp) Successfully created FASTQ files $f1 and ${f2}."
                ((successes++))
            else
                echo "$(time_stamp) Error: $f1 and/or $f2 could not be downloaded." >&2
                ((failures++))
            fi
        else
            f1="${fastq_original_prefix}.fastq"
            if [ -f "$f1" ]; then
                #gzip "$f1"
                write_log "$(time_stamp) Successfully created FASTQ file ${f1}."
                ((successes++))
            else
                write_log "$(time_stamp) Error: $f1 could not be created by fastq-dump."
                ((failures++))
            fi
        fi
        echo >> "$log"
        sleep "$wait_time"
    done
fi

write_log "$(time_stamp) Finished all $accession_count tasks. Downloaded $successes readsets and failed to download $failures readsets."
exit 0
