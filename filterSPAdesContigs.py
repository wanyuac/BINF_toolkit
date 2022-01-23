#!/usr/bin/env python
"""
Filters contigs/scaffolds in SPAdes output FASTA files for a minimum length (bp) and a range of
read depths (min <= d <= max).

Outputs: (1) a filtered FASTA file to stdout, (2) a summary of the filtering process to stderr.
Note that sequence headers in the output FASTA file differ from the original format for the convenience of
subsequent analyses:
	Sequence headers in the input file: NODE_[n]_length_[L]_cov_[C]
	Sequence headers in the output file: NODE_[N] len=[L],cov=[C]
Columns in the output from stderr:
	Input file name, number of contigs passed the filters, number of contigs failed the filters, names of contigs failed the filters

Example command:
	python filterSPAdesContigs.py --input input.fna --min_len 200 --min_d 1 --max_d 100 1>filtered.fna 2>filter.log

Dependencies: Biopython, Python v3

Copyright (C) 2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 23 Jan 2022; the latest update: 23 Jan 2022.
"""
from argparse import ArgumentParser
import os
import sys
from Bio import SeqIO
from collections import namedtuple

def parse_arguments():
	parser = ArgumentParser(description = "Read options and arguments")
	parser.add_argument('--input', '-i', dest = 'input', type = str, required = True, help = "Input FASTA file from SPAdes")
	parser.add_argument('--min_len', '-l', dest = 'min_len', type = int, required = False, default = 1, help = "Minimum contig length [default: 1 bp (no filter)]")
	parser.add_argument('--min_d', '-d0', dest = 'min_d', type = float, required = False, default = 0, help = "Minimum read depth per contig [default: 0 (no filter)]")
	parser.add_argument('--max_d', '-d1', dest = 'max_d', type = float, required = False, default = 0, help = "Maximum read depth per contig [default: 0 (no filter)]")
	return parser.parse_args()

def parse_seq_header(h):
	"""
	Parse sequence headers in SPAdes's output FASTA files
	Format of the headers: NODE_[n]_length_[L]_cov_[C].
	"""
	Contig = namedtuple('Contig', ['name', 'len', 'cov'])
	fields = h.split('_')
	return Contig(name = '_'.join(fields[0 : 2]), len = int(fields[3]), cov = float(fields[5]))

def main():
	args = parse_arguments()
	min_len = args.min_len
	min_d = args.min_d
	max_d = args.max_d
	filter_len = min_len > 1
	filter_min_d = min_d > 0
	filter_max_d = max_d > 0 and min_d < max_d
	fasta = os.path.basename(args.input)
	if not os.path.exists(args.input):
		print(f"Error: input file {fasta} does not exist.", file = sys.stderr)
		sys.exit(1)
	n_pass = 0
	n_fail = 0
	names_fail = []
	for contig in SeqIO.parse(args.input, 'fasta'):
		c = parse_seq_header(contig.id)
		keep = True
		if filter_len:
			keep = keep and c.len >= min_len
		if filter_min_d:
			keep = keep and c.cov >= min_d
		if filter_max_d:
			keep = keep and c.cov <= max_d
		if keep:
			contig.id = c.name
			contig.description = f'len={c.len},cov={c.cov}'
			SeqIO.write(contig, sys.stdout, 'fasta')
			n_pass += 1
		else:
			names_fail.append(contig.id)
			n_fail += 1
	if n_fail > 0:
		ns = ','.join(names_fail)
	else:
		ns = ''
	print(f'{fasta}\t{n_pass}\t{n_fail}\t{ns}', file = sys.stderr)
	
if __name__ == '__main__':
	main()