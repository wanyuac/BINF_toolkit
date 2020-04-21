#!/usr/bin/env python

"""
This script takes a list of NCBI accession numbers (one for each line) from the STDIN and downloads corresponding entries (either GenBank files or FASTA files) under the target directory.

Usage: python download_NCBI_records.py --records "file:objects.txt" --with_prefix --format fasta --email xxx@xxx.com --ext fna --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "NC_0001,NC_0002" --format genbank --email xxx@xxx.com --ext gbk --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "NC_0001,NC_0002" --format genbank --email xxx@xxx.com --prefix K12 --ext gbk --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "file:objects.tsv" --with_prefix --format fasta --email xxx@xxx.com --ext fna --outdir ./ref --skip > download.log
	   Type python download_NCBI_records.py -h or --help for help information.

Important options and arguments:
	--records or -r: Can be either a file (must contain a suffix of ".txt") listing targets to be downloaded, or a string of accession IDs separated by commas (no space is allowed).
	--with_prefix: A logical option specifying that the record file is a tab-delimited file of two columns (without a header line) for accession numbers and prefixes.
	--no_accession: Set this flag to not attach an NCBI accession number after the genome name in each file name. Only applicable when --prefix != None. This option may cause overwriting output files when multiple NCBI accessions share the same prefix.
	--format or -f: The format of files to be downloaded. This option is not used when --db = assembly.
	--db: Customised specification of an NCBI database to retrieve records from
	--ext or -x: The file extension, can be "fasta" (default), "fna", "gb", or "gbk". No dot preceding the extension is needed.
	--outdir or -o: Output directory, no backslash at the end.

An example of the input list: seq_list.txt. Note that accession IDs may not include version numbers, such as ".1".
	HG326223.1\n
	CP011642\n
	
The input file may be composed of two columns: accession number, prefix (genome name), sepearated by a tab character. For instance,
    HG326223.1\tDb11\n
    CP011642\tCAV1492\n

References:
	1. This script is inspired by Mark Schultz's (dr.mark.schultz@gmail.com, GitHub: schultzm) script "downloadGenbankByAccessions.py" stored under the master branch of github.com/katholt/holtlab.
	2. Forum post: www.biostars.org/p/63506/
	3. Damien Farrell's blog post: Retrieving genome assemblies via Entrez with Python (dmnfarrell.github.io/bioinformatics/assemblies-genbank-python)
	
Author: Yu Wan (wanyuac@126.com, https://github.com/wanyuac)
First publication: 27 June 2015 - 14 July 2015; the latest modification: 21 April 2020
Previous name: download_gbk.py
Python version 2 and 3 compatible
Licence: GNU GPL 2.1
"""

from __future__ import print_function
import os
import sys
import time
from Bio import Entrez
from ftplib import FTP
from collections import namedtuple
from argparse import ArgumentParser


def parse_arguments():
	parser = ArgumentParser(description = "Read options and arguments")
	parser.add_argument("--records", "-r", dest = "records", type = str, required = True, \
						help = "Items you want to fetch from the NCBI database")
	parser.add_argument("--with_prefix", "-w", dest = "with_prefix", action = "store_true", required = False, \
						help = "Set when the accession file contains two columns for accessions and prefixes, respectively")
	parser.add_argument("--db", "-d", dest = "db", type = str, default = "nucleotide", required = False, \
						help = "NCBI database to be retrieved from. Options: nucleotide (default), assembly")
	parser.add_argument("--format", "-f", dest = "format", type = str, default = "fasta", required = False, \
						help = "Format: fasta(default)/genbank")
	parser.add_argument("--refseq", "-q", dest = "refseq", action = "store_true", required = False, \
						help = "Set it to specify the RefSeq database for downloading assemblies")
	parser.add_argument("--email", "-e", dest = "email", type = str, required = True, \
						help = "User email address")
	parser.add_argument("--prefix", "-p", dest="prefix", type = str, default = None, required = False, \
						help = "Common prefix adding to all files")
	parser.add_argument("--no_accession", "-n", dest = "no_accession", action = "store_true", required = False, \
						help = "Set this flag to not attach an NCBI accession number after the genome name in each file name. Only applicable when --prefix != None.")
	parser.add_argument("--ext", "-x", dest = "ext", type = str, default = "fasta", required = False, \
						help = "File extension: fasta (default), fna, gb, gbk, fna.gz")
	parser.add_argument("--outdir", "-o", dest = "outdir", type = str, default = ".", required = False, \
						help = "Destination directory, no backslash at the end")
	parser.add_argument("--skip", "-sk", dest = "skip", action = "store_true", required = False, \
						help = "Set to skip downloaded files")
	parser.add_argument("--ftp", "-t", dest = "ftp", type = str, default = "ftp.ncbi.nlm.nih.gov", required = False, \
						help = "Address of the NCBI FTP site from which assemblies are downloaded. Default: ftp.ncbi.nlm.nih.gov.")
	return parser.parse_args()

		
def main():
	args = parse_arguments()
	Entrez.email = args.email
	
	# Read input and set up output file names
	check_output_dir(args.outdir)
	accessions = extract_accessions(targets = args.records, with_prefix = args.with_prefix)
	new_files = create_output_filenames(accessions = accessions, with_prefix = args.with_prefix, \
										no_accession = args.no_accession, outdir = args.outdir, \
										out_prefix = args.prefix, extension = "." + args.ext)
	
	# Iteratively download files
	if args.db == "nucleotide":
		"""
		Download nucleotide records as FASTA files
		"""
		if args.format == "fasta":
			download_records(new_files = new_files, skip_existing = args.skip, record_type = "fasta", \
							 outdir = args.outdir)
		else:
			"""
			Download nucleotide records as GenBank files.
			Do not use "gb" for rettype, as it only includes contig locations if the entry is built
			from contigs.
			"""
			download_records(new_files = new_files, skip_existing = args.skip, record_type = "gbwithparts", \
							 outdir = args.outdir)
	elif args.db == "assembly":
		"""
		Download assemblies as FASTA files
		"""
		download_assemblies(new_files = new_files, skip_existing = args.skip, outdir = args.outdir, \
							use_refseq = args.refseq, site = args.ftp)
	else:
		print("Error: only databases 'nucleotide' and 'assembly' are supported by far. No download task will be launched.")

	return


def extract_accessions(targets, with_prefix):
	"""
	This function deals with two kinds of input: a string of accession numbers or a file.
	"""
	if targets.startswith("file:"):  # treat "targets" as a file name
		with open(targets[5 : ], "r") as f:
			lines = f.read().splitlines()  # read each line into a component of a list and drop newline characters
			
			# parse every line when there are two columns separated by a tab in each line
			if with_prefix:
				accessions = {}
				for line in lines:
					fields = line.split("\t")  # "accession\tprefix"
					accessions[fields[0]] = fields[1]  # return a dictionary
			else:
				accessions = lines  # a list
	else:  # treat "targets" as a series of accession numbers separated by commas.
		accessions = targets.split(",")
		
	return accessions


def create_output_filenames(accessions, with_prefix, no_accession, outdir, out_prefix, extension):
	"""
	Set up file names for download
	"""
	new_files = {}
	if with_prefix:  # when "accessions" is a dictionary
		for entry, prefix in accessions.items():
			if no_accession:
				new_files[entry] = os.path.join(outdir, prefix + extension)
			else:
				new_files[entry] = os.path.join(outdir, prefix + "__" + entry + extension)
	else:
		for entry in accessions:  # when "accessions" is a list
			if out_prefix != None:
				if no_accession:
					new_files[entry] = os.path.join(outdir, out_prefix + extension)
				else:
					new_files[entry] = os.path.join(outdir, out_prefix + "__" + entry + extension)
			else:
				new_files[entry] = os.path.join(outdir, entry + extension)
	
	return new_files


def download_records(new_files, skip_existing, record_type, outdir):
	"""
	Download from the nucleotide database and save files in FASTA or GenBank format
	"""
	print("Start to download records from the NCBI Nucleotide database.")
	n = 0  # the counter for downloaded files
	for entry, new_file in new_files.items():
		if os.path.exists(new_file) and skip_existing:
			print(new_file + " already exists, skipped.")
			continue  # go to the next entry
		try:
			handle = Entrez.efetch(db = "nucleotide", id = entry, rettype = record_type, retmode = "text")
			with open(new_file, "w") as output_file:
				print("Downloading " + entry + " to " + new_file)
				output_file.write(handle.read())  # read and write this entry into a new file named after its accession number
			n += 1
			handle.close()
		except:
			print("The record "+ entry + " is not found.")
			continue
		time.sleep(1) # Pause for one second to obviate submitting too many concurrent requests to NCBI
	print("Done. Altogether %d files were downloaded and stored in %s successfully." % (n, outdir))
	
	return


def download_assemblies(new_files, skip_existing, outdir, use_refseq, site):
	"""
	Download nucleotide sequences from the NCBI Assembly database
	"""
	print("Start to download records from the NCBI Assembly database.")
	urls = get_urls(new_files, use_refseq, skip_existing)  # urls is a named tuple with three fields
	
	# Download files through FTP
	print("Connecting to site " + site + ".")
	try:
		ftp = FTP(site, timeout = 30)  # In general, urlsplit(URL).netloc returns the site address. Paramenter "timeout" is mandatory.
		ftp.login()  # '230 Anonymous access granted, restrictions apply'
		print("Successfully logged in.")
	except:
		sys.exit("Error: cannot log in to the site.")
	
	n = 0
	prefix_len = len("ftp://" + site)  # For example, len("ftp://ftp.ncbi.nlm.nih.gov") = 26
	
	for assembly in urls:
		try:
			with open(assembly.local, "wb") as f: # Create a binary file
				"""
				The retrbinary method of an FTP object does not work for the full FTP address, such as
				ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/.../file.fna.gz. It only accepts the path
				/genomes/all/GCA/000/.../file.fna.gz. Hence we need to extract the path from an URL
				using a simple command: assembly.url[prefix_len : ].
				"""
				ftp.retrbinary("RETR " + assembly.url[prefix_len : ], f.write)
			print("Saved %s as %s." % (assembly.url, assembly.local))
			n += 1
			time.sleep(1)
		except:
			print("Warning: remote file " + assembly.url + " is not accessible. Skip.")
			os.system("rm " + assembly.local)
	ftp.quit()
	print("Done. Altogether %d files were downloaded and stored in %s successfully." % (n, outdir))
	
	return


def get_assembly_summary(seq_id):
	"""
	Retrieve details of an assembly under a given ID (the parameter 'id')
	This is a subordinate function of download_assemblies.
	"""
	handle = Entrez.esummary(db = "assembly", id = seq_id, report = "full")
	record = Entrez.read(handle)
	
	return record


def get_urls(new_files, use_refseq, skip_existing):
	"""
	Retrieve FTP addresses of assembly files on NCBI's server.
	"""
	Assembly = namedtuple("Assembly", ["accession", "url", "local"])
	urls = []
	for entry, new_file in new_files.items():
		if os.path.exists(new_file) and skip_existing:
			print(new_file + " already exists, skipped.")
			continue  # go to the next entry
		try:
			handle = Entrez.esearch(db = "assembly", term = entry, retmax = "1")  # There is only one record per accession number.
			record = Entrez.read(handle)
			seq_id = record["IdList"][0]  # Convert this single-element list into a string
			summary = get_assembly_summary(seq_id)
			if use_refseq:
				url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
			else:
				url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
			
			if url == "":
				print("Warning: URL of %s (sequence ID: %s) is not found." % (entry, seq_id))
				continue
			else:
				urls.append(Assembly(accession = entry, \
									 url = os.path.join(url, os.path.basename(url) + "_genomic.fna.gz"), \
									 local = new_file))
		except:
			print("The record "+ entry + " is not found.")
			continue
		time.sleep(1)
		
	return urls


def check_output_dir(outdir):
	"""
	Prepare the output directory
	"""
	if outdir != ".":
		if not os.path.exists(outdir):
			os.system("mkdir " + outdir)
		else:
			print("Output directory " + outdir + " exists.")
	else:
		print("Skipped checking the output directory as it is the current working directory.")
	
	return


if __name__ == "__main__":
	main()
