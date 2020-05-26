#!/usr/bin/env python

"""
This script extracts nucleotide or protein sequences from a GenBank file in accordance with a list of
(locus_tag/gene, feature type) tuples.

Required modules: Python 3, BioPython, argparse, csv.

Usage:
    python extractSeqFromGBK.py --targets [target file] --gbk [GenBank file(s)] > [output file name]
    python extractSeqFromGBK.py --targets [target file] --gbk [GenBank file(s)] --usegene > [output file name]
    python extractSeqFromGBK.py --targets [target file] --gbk [GenBank file(s)] --usegene --aa --extname > [output file name]

Inputs
    1. --gbk: A list of GenBank files. Each filename is a genome name, which will be parsed and put into output sequence headers.
    2. --targets:
        (1) A tab-delimited text file (.tsv) listing selected locus_tags/gene names in the format: [tag]'\t'[feature_type].
        (2) Single-target mode: an ID, preceded by a '^' sign and followed by a feature type with a colon as the delimiter. For example, 'repB,CDS'.
	
    Allowed feature types are: CDS, tRNA, rRNA, tmRNA.
        For example
            SMDB11_RS00910	rRNA
            SMDB11_RS21915	rRNA
            SMDB11_RS00015	CDS
        Or:
            gene1	CDS
            gene2	CDS
            gene3	CDS

Output (to stdout)
    Nucleotide sequences in FASTA format with the header in the format:
        >[sequence ID] [gene name]|[Genome name]|[NCBI nucleotide accession or contig name]|[Coding strand (+/-)]|[Coordinates]|[Coordinate strand (+/-)]|[Locus tag]|[NCBI protein accession/NA]|[Product name]

Example commands
    python extractSeqFromGBK.py --targets loci.tsv --gbk genome1.gb > genome1_genes.fna
    python extractSeqFromGBK.py --targets ^geneA:CDS --gbk *.gbk --usegene > genes.fna
    python extractSeqFromGBK.py --targets genes.tsv --gbk *.gbk --usegene --aa --extname > proteins.faa

Notes
    1. Locus tags are recommended when users want to extract sequences of exact features, since some features in a record
    may share the same gene name. (Check GenBank files before using this script). Nonetheless, users may want to use gene
    names rather than locus tags to include sequences of the same gene name.
    2. Multiple sequences of any target that is shared by different loci will be extracted. For example, two sequences are
    printed for a gene when both sequences share the same gene name.
    3. Note that some features, such as 'source', does not contain a qualifier "locus_tag". A KeyError will arise if call
    qualifiers["locus_tag"] for those features. Moreover, the 'gene' feature, although it shares the same locus_tag with
    its CDS, it does not contain a nucleotide sequence and hence should not be used as a legit type of target features.
    4. Coordinate strand (+/-) always equals '+' when the GenBank file is created by Prokka or NCBI's Prokaryotic Genome
    Annotation Pipeline (PGAP). This scripts presumes that this is the case.

Explanation of warning(s)
    An "IndexError: list index out of range" will arise if the tag list uses Unicode codes.
	
Copyright (C) Yu Wan 2020 <wanyuac@126.com>
Publication: 19 June 2015; latest update: 26 May 2020
Licence: GNU General Public License v3.0
Previous filename: get_gene_seq.py

References
    Mark Schultz, https://github.com/schultzm/parseGenbank_extractGenes.py
    martineau, http://stackoverflow.com/questions/14734604/python-dictionary-of-lists-from-tab-delimited-file
"""

from Bio import SeqIO
import os
import sys
import csv  # User"s Python script name should not be the module name, otherwise, the former will be loaded and cause an error of no attribute loaded.
from argparse import ArgumentParser


def parse_args():
	parser = ArgumentParser(description= "Extract nucleotide/protein sequences from GenBank files")
	parser.add_argument("--targets", "-t", dest = "targets", type = str, required = True, help = "A tab-delimited file listing (locus_tag/gene name, feature type) tuples, or ^tag:type")
	parser.add_argument("--gbk", "-g", dest = "gbk", nargs = "+", type = str, required = True, help = "One or multiple GenBank files")
	parser.add_argument("--usegene", "-u", dest = "usegene", action = "store_true", required = False, help = "A flag enabling the use of gene names rather than locus tags for feature match")
	parser.add_argument("--aa", "-a", dest = "aa", action = "store_true", required = False, help = "Set to print amino acid sequences instead of nucleotide sequences")
	parser.add_argument("--extname", "-x", dest = "extname", action = "store_true", required = False, help = "Set to attach genome names to sequence names, making an extended sequence name")
	
	return parser.parse_args()

		
def main():
	args = parse_args()

	# Read targets
	targets_def = args.targets
	if targets_def.startswith("^"):
		targets_def = targets_def[1 : ]
		tag, feature_type = targets_def.split(":")
		tags = {tag : feature_type}
	else:
		tags = read_table(args.targets)  # read the table and store it as a dictionary

	# Validity check of tags
	if len(tags) == 0:  # if the dictionary 'tags' is empty, then terminate the loop
		print("Warning: no target is read.")
		sys.exit(0)

	search_key = "gene" if args.usegene else "locus_tag"

	for gbk in args.gbk:
		print("Processing %s" % gbk, file = sys.stderr)
		process_gbk(gbk = gbk, tags = tags, search_key = search_key, usegene = args.usegene, get_protein = args.aa,\
			att_name = args.extname, tag_num = len(tags))
	
	return


def process_gbk(gbk, tags, search_key, usegene, get_protein, att_name, tag_num):
	"""
	This function processes a single GenBank file.
	"""
	target_feature_types = set(tags.values())  # Creates a set of feature types (CDS, tRNA, tmRNA, rRNA, etc.) from the dictionary of targets
	if "gene" in target_feature_types:
		print("Error: 'gene' is not a legit feature type.")
		sys.exit(0)
	targets = list(tags.keys())  # Names of targets. For instance, a list of gene names or locus tags.
	g = os.path.splitext(os.path.basename(gbk))[0]  # Remove path and filename extension from the path of the input GenBank file.
	loci_found = 0  # Number of target loci encountered in this GenBank file. This variable is useful when usegene = False.
	continue_search = True

	for contig in SeqIO.parse(gbk, "genbank"):
		"""
		Do not use list(SeqIO.parse(gbk, "genbank")) in order to save memory.
		Object 'contig' belongs to class SeqRecord and corresponds to a LOCUS feature in the GenBank file.
		A GenBank file may be comprised of multiple contigs. The following loop goes through every feature of the contig.
		"""
		if continue_search:  # Go through features of the current contig.
			for f in contig.features:
				if f.type in target_feature_types:  # Skipping unwanted feature types saves time.
					f_qualifier_keys = list(f.qualifiers.keys())
					if search_key in f_qualifier_keys:  # type(f.qualifiers): collections.OrderedDict
						tag_name = f.qualifiers[search_key][0]  # Equals gene name when search_key is 'gene' or locus tag when search_key is 'locus_tag'.
						if tag_name in targets and f.type == tags[tag_name]: # If the true feature type matches the anticipated type, then it is a true discovery.
							strand = "+" if f.strand == 1 else "-"
							start = int(f.location.start) + 1  # An alias for f.location.nofuzzy_start. Conventional start position is 1 bp greater than the Python-style coordinate.
							end = int(f.location.end)  # An alias for f.location.nofuzzy_end

							# Get the sequence
							if get_protein and f.type == "CDS":
								if "translation" in f_qualifier_keys:
									seq = f.qualifiers["translation"][0]  # Type: str
								else:  # It happens when the CDS is a pseudo gene.
									print("Warning: CDS of feature %s in %s does not have a translated sequence." % (tag_name, gbk),\
										file = sys.stderr)
									continue  # Skip the current feature and move to the next one.
							else:
								seq = str(f.extract(contig.seq))

							# Determine the output sequence ID
							if att_name:
								seq_id = "%s.%s" % (tag_name, g)
							else:
								seq_id = tag_name
							
							# Determine the gene name where available
							if "gene" in f_qualifier_keys:
								gene_name = f.qualifiers["gene"][0]
							else:
								gene_name = "NA"

							# Determine protein accession number when available
							if f.type == "CDS" and "protein_id" in f_qualifier_keys:
								protein_accession = f.qualifiers["protein_id"][0]
							else:
								protein_accession = "NA"

							# Get the locus tag
							if usegene and "locus_tag" in f_qualifier_keys:
								locus_tag = f.qualifiers['locus_tag'][0]
							else:
								locus_tag = tag_name

							# Get the product name
							if "product" in f_qualifier_keys:
								product = f.qualifiers["product"][0]
							else:
								product = "NA"
							
							# Print the target sequence
							print(">%s %s|%s|%s|%s|%i-%i|+|%s|%s|%s" % (seq_id, gene_name, g, contig.id, strand, start, end, locus_tag,\
								protein_accession, product))  # print the header
							print(seq)  # extract nucleotide sequence of this feature

							"""
							In order to save time, the for loop is terminated when locus tags (which are unique in every GenBank file)
							are used as search keys and all target locus tags have been found.
							"""
							if not usegene:
								loci_found += 1
								if loci_found == tag_num:  # Do not need to do further search.
									continue_search = False  # To break the outer for loop
									break  # Terminate the current for loop
		else:
			break  # This termination happens when usegene = False and all target locus tags have been found.

	return


def read_table(f):
	"""
	This function reads a tab-delimited file and saves it as a dictionary through using the first column as keys.
	"""
	d = {}  # create an empty dictionary
	with open(f, "r") as csv_file:  # a wrapper for reading a file instead of using open() and f.close()
		csv_reader = csv.reader(csv_file, delimiter = "\t")
		try:
			for row in csv_reader:
				d[row[0]] = row[1]  # use the first column as keys and the second column as values
		except:
			print("Error: cannot read tags values. Your tag file should use ASCII or utf-8 encoding system.")
			sys.exit(1)
	return d

					
if __name__ == "__main__":
	main()
