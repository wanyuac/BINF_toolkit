'''
This script generates a list of file names based on a list of strings. It is useful if you want to generate a list of file names for read sets from a list of bacterial strain names.

Usage: python filename_generator.py -i <input file> -o <output file> -p <prefix> -s <suffix> -f <from> -l <to> -pe

Input: a list of filenames
	Example: (inlist.txt)
		sample1__genes__results.txt
		sample2__genes__results.txt
Command: python filename_generator.py -i inlist.txt -o outlist.txt -p /reads/ -s .fastq.gz -f 0 -l 7 -pe
Output: a list of new file names generated on the basis of strings in inlist.txt
	Example: (outlist.txt)
		/reads/sample1_1.fastq.gz
		/reads/sample1_2.fastq.gz
		/reads/sample2_1.fastq.gz
		/reads/sample2_2.fastq.gz
		
Author: Yu Wan (wanyuac@gmail.com, GitHub: https://github.com/wanyuac)
First edition: 6 July 2015
Last edition: 5 Nov 2015

License: GNU GPL 2.1
'''

from argparse import ArgumentParser

def parse_args():
	# Read arguments from the command line
	parser = ArgumentParser(description='Regenerate filenames.')
	# Inputs
	parser.add_argument('-i', type = str, required = True, help = 'File name of the input list')
	parser.add_argument('-o', type = str, required = True, help = 'File name of the output list')
	parser.add_argument('-p', type = str, required = False, help = 'The prefix added to the base for new filenames')
	parser.add_argument('-s', type = str, required = False, default = '.fastq.gz', help = 'The suffix added to the base for new filenames')
	parser.add_argument('-f', type = int, required = True, help = 'From which character of the base')
	parser.add_argument('-l', type = int, required = True, help = 'How many characters of the base should be used; -1: use the whole base')
	parser.add_argument('-pe', required = False, action='store_true', help = 'Whether read sets are paired-end')
	return parser.parse_args()

def main():
	args = parse_args()
	with open(args.i, 'rU') as in_f:
		bases = in_f.read().splitlines()
	out_f = open(args.o, 'w')
	
	if args.l > -1:  # only use part of the base for constructing a new file name
		for i in range(0, len(bases)):
			bases[i] = bases[i][args.f : args.l]
		
	for item in bases:
		if args.pe:  # if input files are related to paired-ended libraries
			for i in range(1, 3):
				filename = '{prefix}{base}_{index}{suffix}\n'.format(prefix = args.p, base = item, index = i, suffix = args.s)
				out_f.write(filename)
		else:
			filename = args.p + item[args.f : args.l] + args.s + '\n'
			out_f.write(filename)
	out_f.close()
	print 'All filenames were generated from bases.'

if __name__ == '__main__':
	main()