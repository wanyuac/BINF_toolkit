"""
This script takes a list of NCBI accession numbers (one for each line) from the STDIN and downloads corresponding entries (either GenBank files or FASTA files) under the target directory.

Usage: python download_NCBI_records.py --records "file:objects.txt" --format fasta --email xxx@xxx.com --suffix fna --outdir ./ref --skip > download.log
	   python download_NCBI_records.py --records "NC_0001,NC_0002" --format genbank --email xxx@xxx.com --suffix gbk --outdir ./ref --skip > download.log
	   Type python download_NCBI_records.py -h or --help for help information.

Notes about options and option arguments:
	--records: can be either a file (must contain a suffix of ".txt") listing targets to be downloaded, or a string of accession IDs separated by commas (no space is allowed).
	--format or -f: the format of files to be downloaded
	--suffix or -s: the file extension, can be "fasta" (default), "fna", "gb", or "gbk". No dot preceding the extension is needed.
	--outdir or -o: output directory, no backslash at the end.
An example of the input list: seq_list.txt. Note that accession IDs may not include version numbers, such as ".1".
	HG326223.1\n
	CP011642\n

References:
	1. This script is inspired by Mark Schultz's (dr.mark.schultz@gmail.com, GitHub: schultzm) script "downloadGenbankByAccessions.py" stored under the master branch of https://github.com/katholt/holtlab.
	2. Forum post: www.biostars.org/p/63506/
	
Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Development history
	27 June 2015 - 14 July 2015, 4 November 2015, 6/2/2016
	Previous name: download_gbk.py
Python version: 2.7.5
Licence: GNU GPL 2.1
"""

from Bio import Entrez
from argparse import ArgumentParser
import time, os

def parse_arguments():
	parser = ArgumentParser(description="Read options and arguments")
	parser.add_argument("--records", "-r", dest="records", type=str, required=True, help="Items you want to fetch from NCBI database")
	parser.add_argument("--format", "-f", dest="format", type=str, default="fasta", required=True, help="Format: fasta(default)/genbank")
	parser.add_argument("--email", "-e", dest="email", type=str, required=True, help="User email address")
	parser.add_argument("--suffix", "-s", dest="suffix", type=str, default="fasta", required=False, help="File extension: fasta (default), fna, gb, or gbk")
	parser.add_argument("--outdir", "-o", dest="outdir", type=str, default=".", required=False, help="Destination directory, no backslash at the end")
	parser.add_argument("--skip", "-sk", dest="skip", action="store_true", required=False, help="Set to skip downloaded files")
	return parser.parse_args()

def fetch_records(extension, targets, format, outdir, skip):
	n = 0  # the counter for downloaded files
	
	# deal with two models of input
	if targets.startswith("file:"):  # treat "targets" as a file name
		with open(targets[5 : ], "rU") as f:
			accessions = f.read().splitlines()  # read each line into a component of a list and drop newline characters
	else:  # treat "targets" as a series of accession numbers separated by commas.
		accessions = targets.split(",")
		
	if not os.path.exists(outdir):
		os.system("mkdir " + outdir)
		
	for entry in accessions:
		new_file = os.path.join(outdir, entry + extension)
		if os.path.exists(new_file) and skip:
			print new_file + " is existent, skip."
			continue  # go to the next entry
		try:
			if format == "fasta":
				handle = Entrez.efetch(db = "nucleotide", id = entry, rettype = "fasta", retmode = "text")
			else:
				handle = Entrez.efetch(db = "nucleotide", id = entry, rettype = "gbwithparts", retmode = "text")  # not "gb", which only includes contig locations if the entry is built from contigs
			with open(new_file, "w") as output_file:
				print "Downloading " + entry + " to " + new_file
				output_file.write(handle.read())  # read and write this entry into a new file named after its accession number
			n += 1
			handle.close()
		except:
			print "The record "+ entry + " is not found."
			continue
		time.sleep(1) # should not send too many requests per second to NCBI
	return n
		
def main():
	args = parse_arguments()
	Entrez.email = args.email
	extension = "." + args.suffix  # set the file extension
	count = fetch_records(extension, args.records, args.format, args.outdir, args.skip)
	print "Done. Altogether %d files were downloaded and stored in %s successfully." % (count, args.outdir)

if __name__ == "__main__":
	main()
