'''
This script converts a GenBank file (.gbk or .gb) from Stdin into a Sequin feature table (.tbl), which is an input file of tbl2asn used for creating an ASN.1 file (.sqn).

Package requirement: BioPython and argparse

Usage: python gbk2tbl.py --mincontigsize 200 --prefix <prefix> --modifiers <modifier file> < annotation.gbk 2> stderr.txt

Inputs:
	A GenBank file, which ought to be passed to the script through stdin.
	A modifier file: a plain text file containing modifiers for every FASTA definition line.
		All modifiers must be written in a single line and are separated by a single space character.
		No space should be placed besides the '=' sign. Check http://www.ncbi.nlm.nih.gov/Sequin/modifiers.html for choosing a proper format for modifiers.
		For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]
		This line will be copied and printed along with the record name as the definition line of every contig sequence.

Outputs
	<prefix>.tbl: the Sequin feature table
	<prefix>.fsa: the corresponding fasta file
	These files are inputs for tbl2asn which generates ASN.1 files (.sqn).

Arguments
	--mincontigsize: the minimum contig size, default = 200 in accordance with NCBI's regulation
	--prefix: the prefix of output filenames, default = 'seq'
<<<<<<< HEAD
	--modifiers: Modifiers for every FASTA definition line. All modifiers must be written in a single line and are separated by a single space character.
	  No space should be placed besides the '=' sign. Check http://www.ncbi.nlm.nih.gov/Sequin/modifiers.html for choosing a proper format for modifiers.

Development notes
	This script is derived from the one developed by SEQanswers users nickloman (https://gist.github.com/nickloman/2660685/genbank_to_tbl.py) and ErinL who modified nickloman's script and put it
	on the forum post at http://seqanswers.com/forums/showthread.php?t=19975.
=======
	--modifiers: the filename of the modifier file, default = 'modifiers.txt'
	  
Edition notes
	This script is derived from the one developed by SEQanswers users nickloman (https://gist.github.com/nickloman/2660685/genbank_to_tbl.py) and ErinL who modified nickloman's script and put it
	on the forum post (http://seqanswers.com/forums/showthread.php?t=19975).
	Edition history: 20 June 2015 - 3 July 2015 by Yu Wan (wanyuac@gmail.com)
>>>>>>> origin/master

Author of this version: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Edition history: 20 June and 11 July 2015

Licence: GNU GPL 2.1

Notes for the FASTA header modifiers
	[topology=?]: the molecular topology (circular/linear) of the sequence if this information is not contained in records
		contigs: linear (the default value)
		finished genomes of plasmids and bacterial chromosomes: circular
	An example of the content of the modifier file:
		[organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1]
'''

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

allowed_qualifiers = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag', 'note', 'inference', 'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country']
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
			print >> sys.stderr, 'skipping small contig %s' % (rec.id)
			continue  # start a new 'for' loop
		contig_num += 1
		print rec.name  # print the contig name to STDOUT
		
		# write the fasta file 
		rec.description = modifiers
		SeqIO.write([rec], fasta_fh, 'fasta')  # Prints this contig's sequence to the fasta file. The sequence header will be rec.description.

		# write the feature table
		print >> feature_fh, '>Feature %s' % (rec.name)  # write the first line of this record in the feature table: the LOCUS name
		for f in rec.features:
			# print the coordinates
			if f.strand == 1:
				print >> feature_fh, '%d\t%d\t%s' % (f.location.nofuzzy_start + 1, f.location.nofuzzy_end, f.type)
			else:
				print >> feature_fh, '%d\t%d\t%s' % (f.location.nofuzzy_end, f.location.nofuzzy_start + 1, f.type)

			if (f.type == 'CDS') and ('product' not in f.qualifiers):
				f.qualifiers['product'] = 'hypothetical protein'
			# print qualifiers (keys and values)
			for (key, values) in f.qualifiers.iteritems():
				'''
				Apply the iteritems() method of the dictionary f.qualifiers for (key, values) pairs
				iteritems() is a generator that yields 2-tuples for a dictionary. It saves time and memory but is slower than the items() method.
				'''
				if key not in allowed_qualifiers:
					continue  # start a new 'for' loop of f, skipping the following 'for' statement of v
				for v in values:  # else, write all values under this key (qualifier's name)
					print >> feature_fh, '\t\t\t%s\t%s' % (key, v)
	fasta_fh.close()  # finish the generation of the FASTA file
	feature_fh.close()  # finish the generation of the feature table
	print str(contig_num) + ' records have been converted.'

# call the main function
if __name__ == '__main__':
	main()
