#!/usr/bin/env python

'''
This script adds a sample name at the beginning of each sequence in a FASTA file. For example, the header ">g1 description" becomes
">sample1__g1 description" after running this script.

Author: Yu Wan (wanyuac@gmail.com, github.com/wanyuac)

Example: python add_sample_name_FASTA.py -i filename.txt (or filename.fna) -o output_dir -n

License: GNU GPL 2.0

First edition: Fri 27 Nov 2015
Last edition: Sat 28 Nov 2015
'''

from argparse import ArgumentParser
from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
	parser = ArgumentParser(description="Add a sample name to every header of sequences")
	parser.add_argument("-i", type = str, required = True, help = "A textual list of input files or the file name of a single FASTA file")
	parser.add_argument("-o", type = str, required = False, default = ".", help = "Output directory")
	parser.add_argument("-n", required = False, action="store_true", help = "Whether to extract the sample name from the file name rather than the path?")
	return parser.parse_args()

def main():
	args = parse_args()
	
	# read file names from a list
	if ".txt" in args.i:
		with open(args.i, "rU") as f:
			fasta_files = f.read().splitlines()
	else:
		fasta_files = [args.i]  # If there is just a single FASTA file to be processed.
		
	# read every FASTA file, change all sequence IDs and write into a new file
	for f in fasta_files:
		new_fasta = []
		fields = f.split("/")
		if args.n:
			sample = (fields[-1].split("__"))[0]  # get the sample name from the first part of the file name
		else:
			sample = fields[-3]  # split the path and get the second last field as the sample name
		extension = (fields[-1].split("."))[1]  # get the filename extension: faa, fna or ffn
		records = list(SeqIO.parse(open(f, "rU"), "fasta"))  # records of a single GenBank file
		
		# process each sequence
		for s in records:
			s.id = "__".join([sample, s.id])
			s.description = " ".join(s.description.split(" ")[1 : ])  # remove the first field, which is identical to the sequence ID
			new_fasta.append(SeqRecord(s.seq, id = s.id, name = "", description = s.description))
		
		SeqIO.write(new_fasta, "%s/%s.%s" % (args.o, sample, extension), "fasta")

if __name__ == "__main__":
	main()