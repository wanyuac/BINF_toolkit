#!/usr/bin/python

"""
This script converts a GenBank file (.gbk or .gb) from Stdin into a Sequin feature table (.tbl), which is an input file of tbl2asn used for creating an ASN.1 file (.sqn).

Package requirement: BioPython and argparse

Usage:
	Simple command:
		python gbk2tbl.py --mincontigsize 200 --prefix any_prefix --modifiers modifier_file.txt < annotation.gbk
		cat annotation.gbk | python gbk2tbl.py --mincontigsize 200 --prefix any_prefix --modifiers modifier_file.txt  # integrate gbk2tbl into a pipeline
	Redirecting error messages to a text file (optional):
		python gbk2tbl.py --mincontigsize 200 --prefix any_prefix --modifiers modifier_file.txt < annotation.gbk 2> stderr.txt
		cat annotation.gbk | python gbk2tbl.py --mincontigsize 200 --prefix any_prefix --modifiers modifier_file.txt 2> stderr.txt
	Note that this script reads the GenBank file through the stdin ("< annotation.gbk") and you may want to redirect the stderr to a file via "> stderr.txt" (redirection).
	
Inputs:
	A GenBank file, which ought to be passed to the script through the standard input (stdin).
	A modifier file: a plain text file containing modifiers for every FASTA definition line.
		All FASTA header modifiers must be written in a single line and are separated by a space character. This line will
		be copied and directly printed along with the record name as the definition line of every contig sequence.
		No space should be placed besides the '=' sign. Check http://www.ncbi.nlm.nih.gov/Sequin/modifiers.html for choosing a proper format for modifiers.
		For example, the content of a modifier file can be (no tab character):
			[organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]
		Furthermore, regarding the modifier 'topology':
			[topology=?]: the molecular topology (circular/linear) of the sequence if this information is not contained in records
				For contigs: linear (the default value)
				For finished genomes of plasmids and bacterial chromosomes: circular

Outputs:
	any_prefix.tbl: the Sequin feature table
	any_prefix.fsa: the corresponding fasta file
	These files are inputs for tbl2asn which generates ASN.1 files (*.sqn).

Arguments:
	--mincontigsize: the minimum contig size, default = 200 in accordance with NCBI's regulation
	--prefix: the prefix of output filenames, default = 'seq'
	--modifiers: the filename of the modifier file, default = 'modifiers.txt'
	  
Development notes:
	This script is derived from the one developed by SEQanswers users nickloman (https://gist.github.com/nickloman/2660685/genbank_to_tbl.py) and ErinL who modified nickloman's script and put it
	on the forum post (http://seqanswers.com/forums/showthread.php?t=19975).

Author of this version: Yu Wan (wanyuac@gmail.com, github.com/wanyuac)
Creation: 20 June 2015 - 11 July 2015; the latest edition: 21 October 2019

Dependency: Python versions 2 and 3 compatible.

Licence: GNU GPL 2.1
"""

from __future__ import print_function
import sys
from Bio import SeqIO
from argparse import ArgumentParser

def parse_args():
# Extract arguments from the command line
	parser = ArgumentParser(description= 'Read arguments: species, strain, BioProject, prefix')
	parser.add_argument('--mincontigsize', type = int, required = False, default = 200, help = 'The minimum contig length')
	parser.add_argument('--prefix', type = str, required = False, default = 'seq', help = 'The prefix of output filenames')
	parser.add_argument('--modifiers', type = str, required = True, default = 'modifiers.txt', help = 'The text file containing a single line of FASTA head modifiers')
	return parser.parse_args()

def read_modifiers(file):
# This function only reads the first line of the modifier file. So please ensure that all modifiers are put in the first line.
	with open(file, 'rU') as f:
		s = f.readline()  # only read once
	return s

allowed_qualifiers = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag', 'note', 'inference', \
					  'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country', \
					  'collection_date']  # In GenBank files, the qualifier 'collection-date' is written as 'collection_date'.
'''
These are selected qualifiers because we do not want to see qualifiers such as 'translation', 'transl_table', or 'codon_start' in the feature table.
Qualifiers 'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country' belong to the feature 'source'.
'''

def main():
	args = parse_args()  # read arguments
	contig_num = 0
	fasta_fh = open(args.prefix + '.fsa', 'w')  # the file handle for the fasta file
	feature_fh = open(args.prefix + '.tbl', 'w')  # the file handle for the feature table
	modifiers = read_modifiers(args.modifiers)  # read the modifiers from a text file
	records = list(SeqIO.parse(sys.stdin, 'genbank'))  # read a GenBank file from the standard input and convert it into a list of SeqRecord objects

	for rec in records:  # for every SeqRecord object in the list 'records'
		if len(rec) <= args.mincontigsize:  # filter out small contigs
			print('skipping small contig %s' % (rec.id), file=sys.stderr)
			continue  # start a new 'for' loop
		contig_num += 1
		print(rec.name)  # print the contig name to STDOUT
		
		# write the fasta file 
		rec.description = modifiers
		SeqIO.write([rec], fasta_fh, 'fasta')  # Prints this contig's sequence to the fasta file. The sequence header will be rec.description.

		# write the feature table
		print('>Feature %s' % (rec.name), file = feature_fh)  # write the first line of this record in the feature table: the LOCUS name
		for f in rec.features:
			# print the coordinates
			if f.strand == 1:
				print('%d\t%d\t%s' % (f.location.nofuzzy_start + 1, f.location.nofuzzy_end, f.type), file = feature_fh)
			else:
				print('%d\t%d\t%s' % (f.location.nofuzzy_end, f.location.nofuzzy_start + 1, f.type), file = feature_fh)

			if (f.type == 'CDS') and ('product' not in f.qualifiers):
				f.qualifiers['product'] = 'hypothetical protein'
			# print qualifiers (keys and values)
			for (key, values) in f.qualifiers.items():
				'''
				Apply the iteritems() method of the dictionary f.qualifiers for (key, values) pairs
				iteritems() is a generator that yields 2-tuples for a dictionary. It saves time and memory but is slower than the items() method.
				'''
				if key not in allowed_qualifiers:
					continue  # start a new 'for' loop of f, skipping the following 'for' statement of v
				for v in values:  # else, write all values under this key (qualifier's name)
					print('\t\t\t%s\t%s' % (key, v), file = feature_fh)
	fasta_fh.close()  # finish the generation of the FASTA file
	feature_fh.close()  # finish the generation of the feature table
	print(str(contig_num) + ' records have been converted.')

# call the main function
if __name__ == '__main__':
	main()
