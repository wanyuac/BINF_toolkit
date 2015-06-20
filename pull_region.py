'''
This script extracts a region of nucleotides by positions from a fasta file
Author: Yu Wan (wanyuac@126.com)
Date: 1 June 2015
GitHub: https://github.com/wanyuac/BINF_toolkit
Licence: GNU GENERAL PUBLIC LICENSE Version 2

Previous name: extract_nc_region.py

Arguments
	-i: the path of the input file
	-n: the number of selected contig
	-s: the first nucleotide to be selected
	-e: the last nucelotide to be selected
	-o: the filename of the output
	
Requirements
	Only one region should be selected
	The start and end positions should not spill out
'''

from argparse import (ArgumentParser, FileType)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
	# This function extracts arguments from the command line
	parser = ArgumentParser(description= "Read arguments: input filename, start position, end position, and out filename")
	parser.add_argument("-i", type = str, required = True, help = "Input path")  # append an argument to variable "parser"
	parser.add_argument("-n", type = int, default = 1, help = "The number of selected contig")
	parser.add_argument("-s", type = int, required = True, help = "The start position")
	parser.add_argument("-e", type = int, required = True, help = "The end position")
	parser.add_argument("-o", type = str, default = "selected_region.fasta", help = "The output filename")
	return parser.parse_args()
	
def main():
	args = parse_args()  # read arguments from the command line
	input = args.i
	n = args.n
	start = args.s
	end = args.e
	length = end - start + 1  # the length of selected region
	output = args.o
	i = 1
	L = 0  # the record length
	print "Arguments:\n" + input + "\n" + str(start) + " " + str(n) + " " + str(end) + "\n" + output
	f = open(input, "rU")  # supports universal newlines
	for record in SeqIO.parse(f, "fasta"):
		if i == n:  # if this is the selected contig
			seq = record.seq[start - 1: end]  # gets the selected sequence of this record
			id = str(record.id) + "|" + str(start) + ".." + str(end) + "|" + str(length) + " bp\n" # gets the header of this record
			r = SeqRecord(seq, id, name = "", description = "")
			L = len(seq)
			f_output = open(output, "w")
			SeqIO.write(r, f_output, "fasta")  # saves the selection
			f_output.close()
			print "Length = " + str(L)
		else:
			i = i + 1
			print "No sequence found!"
	f.close()
	
# The main program
if __name__ == "__main__":
    main()