'''

This script converts a GenBank file (.gbk or .gb) into a Sequin feature table (.tbl), which is an input file of tbl2asn used for creating an ASN.1 file (.sqn).

Package requirement: BioPython and argparse

Usage: python gbk2tbl.py --species "<species name>" --strain <strain name> --mincontigsize 200 --topology <linear/circular> --prefix <prefix> < annotation.gbk 2> stderr.txt

Input: a GenBank file

Outputs
	<prefix>.tbl: the Sequin feature table
	<prefix>.fsa: the corresponding fasta file
	These files are inputs for tbl2asn which generates ASN.1 files (.sqn).

Arguments
	--species: the species name, default = "species"
	--strain: the strain name, default = "strain"
	--mincontigsize: the minimum contig size, default = 200 in accordance with NCBI's regulation.
	--topology: the molecular topology (circular/linear) of the sequence if this information is not contained in records
		contigs: linear
		finished genomes of plasmids and bacterial chromosomes: circular
	--prefix: the prefix of output filenames, default = "seq"

Edition notes
	This script is derived from the one developed by SEQanswers users nickloman (https://gist.github.com/nickloman/2660685/genbank_to_tbl.py) and ErinL who modified nickloman's script and put it
	on the forum post (http://seqanswers.com/forums/showthread.php?t=19975).
	Date of last edition: 20 June 2015 by Yu Wan.

Authors ordered by editions: nickloman, ErinL, and Yu Wan (wanyuac@126.com).

Licence: GNU GPL 2.1

'''

import sys
from Bio import SeqIO
from argparse import (ArgumentParser, FileType)

def parse_args():
# Extract arguments from the command line
	parser = ArgumentParser(description= "Read arguments: species, strain, BioProject, prefix")
	parser.add_argument("--species", type = str, required = False, default = "species", help = "The species name")
	parser.add_argument("--strain", type = str, required = False, default = "strain", help = "The strain name")
	parser.add_argument("--mincontigsize", type = int, required = False, default = 200, help = "The minimum contig length")
	parser.add_argument("--topology", type = str, required = False, default = "linear", help = "The molecular topology (linear/circular)")
	parser.add_argument("--prefix", type = str, required = False, default = "seq", help = "The prefix of output filenames")
	return parser.parse_args()
	
allowed_qualifiers = ["locus_tag", "gene", "product", "pseudo", "protein_id", "gene_desc", "old_locus_tag", "note", "inference"]
# These are selected qualifiers because we do not want to see qualifiers such as "translation", "transl_table", or "codon_start" in the feature table.

def main():
	args = parse_args()  # read arguments
	contig_num = 0
	fasta_fh = open(args.prefix + ".fsa", "w")  # the file handle for the fasta file
	feature_fh = open(args.prefix + ".tbl", "w")  # the file handle for the feature table
	records = list(SeqIO.parse(sys.stdin, "genbank"))  # read a GenBank file from the standard input and convert it into a list of SeqRecord objects

	for rec in records:  # for every SeqRecord object in the list "records"
		if len(rec) <= args.mincontigsize:  # filter out small contigs
			print >> sys.stderr, "skipping small contig %s" % (rec.id)
			continue  # start a new "for" loop
		contig_num += 1
		print rec.name
		
		topo = rec.annotations.get("molecule", args.topology)  # returns either "circular" or "linear" if the key "molecule" is present
		'''
		SeqRecord.annotations is a dictionary, hence the dictionary method get(key, default value) applies.
		However, "molecule" is not a key yet, so this method always returns the default value args.topology.
		SeqRecord.annotations may include the key "molecule" in future versions.
		'''
		rec.description = "[organism = %s] [strain = %s] [molecule = DNA] [topology = %s] [tech=wgs] [gcode=11]" % (args.species, args.strain, topo)
		SeqIO.write([rec], fasta_fh, "fasta")  # add the sequence of this contig to the fasta file

		print >> feature_fh, ">Feature %s" % (rec.name)  # write the first line of this record in the feature table: the LOCUS name
		for f in rec.features:
			# print the coordinate
			if f.strand == 1:
				print >> feature_fh, "%d\t%d\t%s" % (f.location.nofuzzy_start + 1, f.location.nofuzzy_end, f.type)
			else:
				print >> feature_fh, "%d\t%d\t%s" % (f.location.nofuzzy_end, f.location.nofuzzy_start + 1, f.type)

			if (f.type == "CDS") and ("product" not in f.qualifiers):
				f.qualifiers["product"] = "hypothetical protein"
			# print qualifiers (keys and values)
			for (key, values) in f.qualifiers.iteritems():
				'''
				Apply the iteritems() method of the dictionary f.qualifiers for (key, values) pairs
				iteritems() is a generator that yields 2-tuples for a dictionary. It saves time and memory but is slower than the items() method.
				'''
				if key not in allowed_qualifiers:
					continue  # start a new "for" loop of f, skipping the following "for" statement of v
				for v in values:  # else, write all values under this key (qualifier's name)
					print >> feature_fh, "\t\t\t%s\t%s" % (key, v)
	print str(contig_num) + " records have been converted."

# call main function
if __name__ == '__main__':
	main()
