#!/bin/bash
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 3/3/2020

# Guidance #########################
display_usage(){
    echo "
    Usage:
        Rename files according to a two-column TSV file: [old filename]\t[new filename]
        Example command: renameFiles.sh names.tsv
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Implementation #########################
while read line; do
    # Split the delimited string into an arrary of two elements.
    # Do not use IFS=$"\t" as it does not work correctly.
    # # https://unix.stackexchange.com/questions/410710/splitting-a-line-into-array-in-bash-with-tab-as-delimiter
    
    IFS=$'\t' read -r -a names <<< "$line"
    echo -e "Change ${names[0]} into ${names[1]}."
    mv ${names[0]} ${names[1]}
done < "$1"  # expect a file name as an input

