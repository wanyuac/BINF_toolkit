#!/bin/bash
# Run this script in conda environment 'ariba'
# Yu Wan (1/7/2021)

cd $1

while read -r g
do
	f="ariba_out/${g}/report.tsv"
	if [ -f "$f" ]
	then	
		cp $f report/in.report.${g}.tsv
	else
		echo "Isolate $g did not have any report generated."
    fi
done < "$2"

ariba summary --no_tree --verbose summary report/in.report.*.tsv
