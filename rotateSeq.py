#!/usr/bin/env python
"""
This script rotates circular sequences and makes them to start from specific positions in output FASTA files. It also circularises every sequence
and save it in a GFA file if a filename is given to argument --gfa/-g. This script is useful for improving complete genome assemblies as well as
read simulation.

Command:
	python rotateSeq.py -i input.fna -t new_starts.tsv -f output.fna -g output.gfa 2> messages.err
	python rotateSeq.py -i input.fna -t new_starts.tsv -f output.fna && gzip output.fna

New start positions are specified in a three-column, tab-delimited, header-free table (parameter '-t'): 'sequence ID'\t'Position'\t'Orientation (+/-)'.
No change is applied to input sequences whose sequence IDs are not listed in this table (e.g., when some sequences are linear or incomplete).

Dependencies: Python 3, BioPython 1.78+.

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 17 June 2021; the latest update: 9 January 2022
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq  # Bio.Alphabet has been removed from BioPython from v1.78. See https://biopython.org/wiki/Alphabet.
from argparse import ArgumentParser
from collections import namedtuple


def parse_argument():
	parser = ArgumentParser(description = "Restart circular sequences from given positions")
	parser.add_argument('-i', '--input', dest = 'i', type = str, required = True, help = "An input FASTA file")
	parser.add_argument('-t', '--table', dest = 't', type = str, required = True, \
	help = "A tab-delimited, header-free table of three columns: sequence ID, the base to be used as the first base of the new sequence, and the orientation (+/-) of the new sequence")
	parser.add_argument('-f', '--fasta', dest = 'f', type = str, required = False, default = 'rotated_seq.fna', help = "Output FASTA file")
	parser.add_argument('-g', '--gfa', dest = 'g', required = False, default = None, help = "Output GFA file, assuming all contigs are circular")
	return parser.parse_args()


def main():
	args = parse_argument()
	prev_seqs = import_seqs(args.i)
	pos_spec = import_positions(args.t)
	fasta_out = open(file = args.f, mode = 'w', encoding = 'utf-8')
	make_gfa = args.g != None
	if make_gfa:
		gfa_out = open(file = args.g, mode = 'w', encoding = 'utf-8')
		gfa_out.write("H\tVN:Z:1.0\n")  # The GFA header
	for i, contig in prev_seqs.items():
		if i in pos_spec.keys():
			p = pos_spec[i]  # 'p' is a namedtuple 'Pos' defined by function import_positions.
			if p.ori == '+':  # Simpler case: start from the chosen base and then go clockwise to create the new sequence.
				if p.base > 0:  # Convert the p-th position into Python's character index. Note that no change will be carried out by Python if p > len(contig.seq).
					s = p.base - 1  # Index of the chosen base (start base) in the sequence (a string).
					print("Restart sequence %s from base %i and go clockwise to create the new sequence." % (i, p.base), file = sys.stderr)
					seq = str(contig.seq)
					contig.seq = Seq(seq[s : ] + seq[ : s])  # Rotate the current sequence; no change when p.base = 0. "generic_dna" is no longer needed from BioPython v1.78.
				else:
					print("Warning: position %i for sequence %s cannot be negative. No change will be applied to this sequence." % (p.base, i), file = sys.stderr)   
			else:
				if p.base > 0:
					s = p.base  # The calculation of the index is: s = (p.base + 1) - 1, where 'p.base + 1' is the new start base required for creating the correct sequence that goes counterclockwise.
					print("Restart sequence %s from base %i and go counterclockwise to create the new sequence." % (i, p.base), file = sys.stderr)
					seq = str(contig.seq)
					contig.seq = Seq(seq[s : ] + seq[ : s]).reverse_complement()  # Rotate and then take the reverse complement (return value: a new Seq object)
				else:
					print("Warning: position %i for sequence %s cannot be negative. No change will be applied to this sequence." % (p.base, i), file = sys.stderr)
		else:
			print("Warning: sequence " + i + " is not found in the position table. No change will be applied to this sequence.", file = sys.stderr)
		SeqIO.write(contig, fasta_out, "fasta")
		if make_gfa:
			gfa_out.write(fasta2gfa(i, contig))
	fasta_out.close()
	if make_gfa:
		gfa_out.close()
	return


def import_seqs(fasta):
	""" Returns a dictionary of SeqIO.SeqRecord objects {seq_id : SeqRecord}. """
	check_file(fasta)
	seqs = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
	return seqs


def import_positions(tsv):
	"""
	Returns a dictionary of namedtuples: {seq_id : Pos(base, ori)}.
		base: he base to be used as the first base of the new sequence;
		ori: to go clockwise (+) or counterclockwise (-) from the chosen base to construct the new sequence.
	"""
	check_file(tsv)
	Pos = namedtuple('Pos', ['base', 'ori'])
	with open(tsv, 'r') as f:
		pos = dict()  # {seq_id : Pos}
		lines = f.read().splitlines()
		for line in lines:
			if line != '':
				i, p, d = line.split('\t')
				pos[i] = Pos(base = int(p), ori = d)
	return(pos)


def fasta2gfa(i, contig):
	""" Makes GFA-formatted lines from an input contig sequence """
	s = str(contig.seq)
	bp = str(len(s))
	lines = '\t'.join(['S', i, s, 'LN:i:' + bp, 'RC:i:' + bp]) + '\n'  # Arbitrarily assign a (relative) fold coverage of one to each sequence.
	lines += '\t'.join(['L', i, '+', i, '+', '0M']) + '\n'
	return lines


def check_file(f):
	if not os.path.exists(f):
		print("Argument error: file " + f + " is not accessible.", file = sys.stderr)
		sys.exit(1)
	return


if __name__ == '__main__':
	main()
