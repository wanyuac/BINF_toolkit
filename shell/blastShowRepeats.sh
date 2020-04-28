#!/usr/bin/bash
# We BLAST a nucleotide sequence against itself to identify repetitive regions. Of course, every region
# matches to itself as well. This script is hence developed to remove such self-matches. It assumes
# the input crunch file follows default columns in the '-fmt 6' output format of BLAST:
#    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#
# Example commands (First, use the 'chmod' command to make this script executable):
#   blastn -query sample.fasta -db ref -task megablast -evalue 0.01 -perc_identity 90 -max_target_seqs 10 -outfmt 6 | blastShowRepeats.sh sample_vs_ref.crunch
#   or: blastShowRepeats.sh sample_vs_ref.crunch > sample_vs_ref_repeats.crunch
#
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 28 Apr 2020

# Parse tab-delimited lines into an array
while IFS=$'\t' read -r -a line
do
    qstart="${line[6]}"
    qend="${line[7]}"
    sstart="${line[8]}"
    send="${line[9]}"
    
    # The following statement is the same as [[ "$qstart" -ne "$sstart" && "$qend" -ne "$send") ]].
    # This if statement ignores hits where qstart = sstart AND qend = send.
    if [ "$qstart" -ne "$sstart" ] && [ "$qend" -ne "$send" ]
    then
        # We use a sub-shell to avoid overriding the current IFS: ( IFS=$'\t'; echo "${line[*]}" ).
        # https://superuser.com/questions/461981/how-do-i-convert-a-bash-array-variable-to-a-string-delimited-with-newlines/462400
        # It is also necessary here even though the IFS has been set to a tab character for read data.
        # Otherwise, the output of echo is space-delimited.
        ( IFS=$'\t'; echo "${line[*]}" )
    fi
done < "${1:-/dev/stdin}"
