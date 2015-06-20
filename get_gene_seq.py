"""
This script extracts gene sequences from a GenBank file, in accordance with a list of (locus_tag, feature type) tuples.
Required module: Bio, argparse, csv

Usage: python get_gene_seq.py --tags <locus_tag file> --gb <GenBank file> > genes.fna

Inputs
	A GenBank file.
	A text file listing selected locus_tags in the following format: locus_tag"\t"feature_type. This file MUST use ASCII codes because the module csv/2.3 does not support Unicode inputs (https://docs.python.org/2/library/csv.html).
	Allowed feature types are: CDS, tRNA, rRNA, tmRNA.
		For example
			SMDB11_RS00910	rRNA
			SMDB11_RS21915	rRNA
			SMDB11_RS00015	CDS
Output
	Nucleotide sequences in FASTA format with the header in the format:
		>feature type|contig name|locus_tag|position|length|product

Warnings
	1. Although it is unlikely in a GenBank file, but please always ensure that there is no duplication of locus_tag"s in the table because this script treats locus_tag"s as keys for retrieving feature types.
	2. An "IndexError: list index out of range" will arise if the tag list uses Unicode codes.
	
Author: Yu Wan (wanyuac@126.com, GitHub: https://github.com/wanyuac)
Date of last edition: 19 June 2015
Licence: GNU GPL 2.0

References
	Mark Schultz, https://github.com/schultzm/parseGenbank_extractGenes.py
	martineau, http://stackoverflow.com/questions/14734604/python-dictionary-of-lists-from-tab-delimited-file
"""

#! /usr/bin/python

from Bio import SeqIO
import csv  # User"s Python script name should not be the module name, otherwise, the former will be loaded and cause an error of no attribute loaded.
from argparse import (ArgumentParser, FileType)

def parse_args():
# Extract arguments from the command line
	parser = ArgumentParser(description= "Read arguments: tags and gb")
	parser.add_argument("--tags", type = str, required = True, help = "A list of (locus_tag, feature type) tuples")
	parser.add_argument("--gb", type = str, required = True, help = "The GenBank file")
	return parser.parse_args()
	
def read_table(file):
# This function reads a tab delimited file and save it as a dictionary, which uses the first column as keys.
	d = {}  # create an empty dictionary
	with open(file, "rU") as csv_file:  # a wrapper for reading a file instead of using open() and f.close()
		csv_reader = csv.reader(csv_file, delimiter="\t")  # only in ASCII codes
		try:
			for row in csv_reader:
				d[row[0]] = row[1]  # use the first column as keys and the second column as values
		except:
			print "Cannot get tags. Your tag file should use ASCII codes."
			traceback.print_exc()
			raise
	return d
	
target_features = ["CDS", "rRNA", "tRNA", "tmRNA"]  # Only these features are wanted.
		
def main():
	args = parse_args()
	tags = read_table(args.tags)  # read the table and store it as a dictionary
	
	for rec in SeqIO.parse(args.gb, "genbank"):
		if len(tags) == 0:  # if the tags is empty, then terminate the loop
			break
		# If there are still some gene sequences to be extracted, then find every (remaining) locus_tag in each SeqFeature.
		for f in rec.features:
			if len(tags) == 0:  # If all targets have been found, then quit the loop.
				break
			if f.type in target_features:
			# This is important because some features, such as "source", does not contain a qualifier "locus_tag". A KeyError will arise if call qualifiers["locus_tag"] for those features.
			# Moreover, by using target_features, I can skip "gene" features which share the same locus_tag with its CDS but do not contain a nucleotide sequence.
				locus_tag = f.qualifiers["locus_tag"][0]
				if locus_tag in tags.keys() and f.type == tags[locus_tag]: # If the true feature type matches the anticipated type, then it is a true discovery.
					seq = str(f.extract(rec.seq))
					print ">%s|%s|%s|%d|%s" % (rec.name, locus_tag, f.location, len(seq), f.qualifiers["product"][0])  # print the header
					print seq  # extract nucleotide sequence of this feature
					del tags[locus_tag]  # Delete this item from the dictionary because we have already found the gene.
							
if __name__ == "__main__":
	main()