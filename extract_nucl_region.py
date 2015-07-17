"""
This script extracts a region of nucleotides by positions from a fasta file
Author: Yu Wan (wanyuac@126.com)
Date: 1 June and 17 July 2015
GitHub: https://github.com/wanyuac/BINF_toolkit
Licence: GNU GENERAL PUBLIC LICENSE Version 2

Previous name: extract_nc_region.py

Arguments
	-i: the path of the input file
	-n: the name of your selected contig
	-f: feature name specified by the user
	-s: the first nucleotide to be selected
	-e: the last nucelotide to be selected
	-o: the filename of the output
	
Requirements
	Only one region should be selected
	The start and end positions should not spill out
"""

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
	# This function extracts arguments from the command line
	parser = ArgumentParser(description="Read arguments: input filename, start position, end position, and out filename")
	parser.add_argument("-i", type=str, required=True, help="Input path")  # append an argument to variable "parser"
	parser.add_argument("-c", type=str, default="", help="Name of your selected contig; the first contig will be chosen if -c is not set.")
	parser.add_argument("-f", type=str, default="feature", help="The feature name")
	parser.add_argument("-s", type=int, required=True, help="Start position")
	parser.add_argument("-e", type=int, required=True, help="End position")
	parser.add_argument("-o", type=str, default="selected_region.fasta", help="File name of the output")
	return parser.parse_args()
	
def write_seq(contig, feature, start, end, output):
	# read and write FASTA files
	seqlen = end - start + 1  # the length of selected region
	contiglen = len(contig.seq)
	if start > contiglen or end > contiglen:  # the genetic coordinates are out bounded
		flag = False
	else:
		seq = contig.seq[start - 1: end]  # gets the selected sequence of this contig
		descr = feature + "|" + str(start) + ".." + str(end) + "|" + str(seqlen) + " bp\n" # gets the header of this contig
		new_rec = SeqRecord(seq=seq, id=contig.id, name=feature, description=descr)  # create a new SeqRecord instance. Note that the contig.name will not be written in a FASTA file (only in GenBank files).
		f = open(output, "w")
		SeqIO.write(new_rec, f, "fasta")  # saves the selection
		f.close()
		flag = True
	return flag
	
def main():
	args = parse_args()  # read arguments from the command line
	f = open(args.i, "rU")  # supports universal newlines
	contigs = list(SeqIO.parse(f, "fasta"))
	found = False
	if args.c == "":
		found = write_seq(contig=contigs[0], feature=args.f, start=args.s, end=args.e, output=args.o)  # read the first contig if -n is not set
	else:
		for contig in contigs:
			if contig.id == args.c:  # if this is the selected contig
				found = write_seq(contig=contig, feature=args.f, start=args.s, end=args.e, output=args.o)
				break
	f.close()
	if found:
		print "The target sequence was extracted."
	else:
		print "No sequence was found."
	
# The main program
if __name__ == "__main__":
    main()
	