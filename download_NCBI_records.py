"""
This script takes a list of NCBI accession numbers (one for each line) from the STDIN and downloads corresponding entries (either GenBank files or FASTA files) under the target directory.

Usage: python download_NCBI_records.py --records "file:objects.txt" --format fasta --email xxx@xxx.com --ext fna --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "NC_0001,NC_0002" --format genbank --email xxx@xxx.com --ext gbk --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "NC_0001,NC_0002" --format genbank --email xxx@xxx.com --prefix K12 --ext gbk --outdir ./ref --skip > download.log
	   Type python download_NCBI_records.py -h or --help for help information.

Important options and arguments:
	--records or -r: can be either a file (must contain a suffix of ".txt") listing targets to be downloaded, or a string of accession IDs separated by commas (no space is allowed).
	--format or -f: the format of files to be downloaded
	--ext or -x: the file extension, can be "fasta" (default), "fna", "gb", or "gbk". No dot preceding the extension is needed.
	--outdir or -o: output directory, no backslash at the end.
An example of the input list: seq_list.txt. Note that accession IDs may not include version numbers, such as ".1".
	HG326223.1\n
	CP011642\n

References:
	1. This script is inspired by Mark Schultz's (dr.mark.schultz@gmail.com, GitHub: schultzm) script "downloadGenbankByAccessions.py" stored under the master branch of https://github.com/katholt/holtlab.
	2. Forum post: www.biostars.org/p/63506/
	
Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Development history
	27 June 2015 - 14 July 2015, 4 November 2015, 6/2/2016, 23/7/2016, 19/11/2017
	Previous name: download_gbk.py
Python version 2 and 3 compatible
Licence: GNU GPL 2.1
"""

from __future__ import print_function
from Bio import Entrez
from argparse import ArgumentParser
import time, os


def parse_arguments():
	parser = ArgumentParser(description="Read options and arguments")
	parser.add_argument("--records", "-r", dest = "records", type = str, required = True, help = "Items you want to fetch from NCBI database")
	parser.add_argument("--format", "-f", dest = "format", type = str, default = "fasta", required = True, help = "Format: fasta(default)/genbank")
	parser.add_argument("--email", "-e", dest = "email", type = str, required = True, help = "User email address")
	parser.add_argument("--prefix", "-p", dest="prefix", type = str, required = False, default = None, help = "Common prefix adding to all files")
	parser.add_argument("--ext", "-x", dest = "ext", type = str, default = "fasta", required = False, help = "File extension: fasta (default), fna, gb, or gbk")
	parser.add_argument("--outdir", "-o", dest = "outdir", type = str, default = ".", required = False, help = "Destination directory, no backslash at the end")
	parser.add_argument("--skip", "-sk", dest = "skip", action = "store_true", required = False, help = "Set to skip downloaded files")
	return parser.parse_args()

		
def main():
	args = parse_arguments()
	
	extension = "." + args.ext  # set the file extension
	Entrez.email = args.email
	targets = args.records
	n = 0  # the counter for downloaded files
	
	# deal with two kinds of input
	if targets.startswith("file:"):  # treat "targets" as a file name
		with open(targets[5 : ], "rU") as f:
			accessions = f.read().splitlines()  # read each line into a component of a list and drop newline characters
	else:  # treat "targets" as a series of accession numbers separated by commas.
		accessions = targets.split(",")
		
	if not os.path.exists(args.outdir):
		os.system("mkdir " + args.outdir)
		
	for entry in accessions:
		if args.prefix != None:
			new_file = os.path.join(args.outdir, args.prefix + "__" + entry + extension)
		else:
			new_file = os.path.join(args.outdir, entry + extension)
			
		if os.path.exists(new_file) and args.skip:
			print(new_file + " already exists, skipped.")
			continue  # go to the next entry
		try:
			if args.format == "fasta":
				handle = Entrez.efetch(db = "nucleotide", id = entry, rettype = "fasta", retmode = "text")
			else:
				handle = Entrez.efetch(db = "nucleotide", id = entry, rettype = "gbwithparts", retmode = "text")  # not "gb" for rettype, as it only includes contig locations if the entry is built from contigs
			with open(new_file, "w") as output_file:
				print("Downloading " + entry + " to " + new_file)
				output_file.write(handle.read())  # read and write this entry into a new file named after its accession number
			n += 1
			handle.close()
		except:
			print("The record "+ entry + " is not found.")
			continue
		time.sleep(1) # Pause for one second to obviate submitting too many concurrent requests to NCBI
	
	print("Done. Altogether %d files were downloaded and stored in %s successfully." % (n, args.outdir))
	

if __name__ == "__main__":
	main()
