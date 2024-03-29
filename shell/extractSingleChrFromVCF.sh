#!/bin/sh

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 28/4/2021; latest update: 28/4/2021

display_useage() {
    echo "
    Extract variants in a specific chromosome from a VCF file.
    Command: extractSingleChrFromVCF.sh [input VCF] [output VCF] [target chromosome name]
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

vcf_in=$1
vcf_out=$2
chr=$3
outdir=$(dirname $vcf_out)
tmpfile=$outdir/tmp.vcf

n=$(grep -n "##contig=<ID=$chr" $vcf_in | awk -F: '{print $1}')  # https://stackoverflow.com/questions/3213748/get-line-number-while-using-grep
head -n $n $vcf_in > $vcf_out

# The following two lines avoid printing a duplicated "##contig=<ID=$chr ..." line.
((n++))
tail -n +$n $vcf_in > $tmpfile

grep '#CHROM' $tmpfile >> $vcf_out
grep "$chr" $tmpfile >> $vcf_out

rm -f $tmpfile
