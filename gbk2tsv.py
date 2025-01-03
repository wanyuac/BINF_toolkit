#!/usr/bin/python

"""
Convert GenBank files to tab-delimited text files (*.tsv). Every GenBank file must contain the locus_tag qualifier and may
contain multiple contigs (LOCUS).

Usage:
python gbk2tsv.py --gbk 1.gbk --outdir . --features "CDS,rRNA,tRNA" --nucl_seq --prot_seq
python gbk2tsv.py --gbk 1.gbk 2.gbk 3.gbk --outdir . --features "CDS,rRNA,tRNA" --nucl_seq --prot_seq
python gbk2tsv.py --gbk $(ls *.gbk) --outdir . --features "CDS,rRNA,tRNA" --nucl_seq --prot_seq

An example showing columns in every output file:
Contig    Locus         Feature    Start    End    Strand    Pseudo    Product             Gene     Nucl_seq    Prot_seq
Contig_1     locus_tag_1   CDS        1        2200   +      N         dehydrogenase I     unknown  ...         ...
Contig_1     locus_tag_2   CDS        2230     3100   -      Y         homoserine kinase   unknown  ...         ...
...

Dependency: BioPython
Python versions 2 and 3 compatible
Copyright 2019 Yu Wan (wanyuac@sina.cn)
Licensed under the Apache License, Version 2.0
First version: 13 Sep 2019 (Happy Mid-Autumn Festival)
Latest update: 3 Oct 2025
"""


from __future__ import print_function
from __future__ import division
import os
import sys
import glob
from Bio import SeqIO, SeqFeature
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description = "Convert GenBank files to tab-delimited text files")
    parser.add_argument("-g", "--gbk", nargs = "+", type = str, required = True, dest = "gbks", default = "", help = "Input GenBank files")
    parser.add_argument("-o", "--outdir", type = str, required = False, dest = "outdir", default = ".", help = "Output directory (no backslash or forward slash)")
    parser.add_argument("-f", "--features", type = str, required = False, dest = "features", default = "CDS,tRNA,rRNA", help = "Comma-separated features to store (default CDS,tRNA,rRNA)")
    parser.add_argument("-n", "--nucl_seq", action = "store_true", required = False, dest = "nucl_seq", help = "Turn on this option to print nucleotide sequences of features")
    parser.add_argument("-p", "--prot_seq", action = "store_true", required = False, dest = "prot_seq", help = "Turn on this option to print protein sequences of CDS")
    
    return parser.parse_args()  # An instance of the class ArgumentParser


def main():
    args = parse_args()
    gbk_list = get_input_filenames(args.gbks)
    if args.outdir and not os.path.exists(args.outdir):  # The first logical condition becomes "false" if args.outdir = "".
        os.makedirs(args.outdir)
    
    if (len(gbk_list) == 0):
        sys.exit("Invalid --gbk argument: no GenBank file is found.")
        
    header = ["Contig", "Locus", "Feature", "Start", "End", "Strand", "Pseudo", "Product", "Gene"]
    if args.nucl_seq:
        header += ["Nucl_seq"]
    if args.prot_seq:
        header += ["Prot_seq"]
    
    features = args.features.split(",")  # Features of interest
    if len(features) == 0:
        sys.exit("Invalid --features argument: there is no feature to be extracted.")
    
    for gbk in gbk_list:
        tsv_name = os.path.join(args.outdir, os.path.splitext(os.path.basename(gbk))[0] + ".tsv")  # Define the current output filename: pwd/1.gbk -> pwd/1.tsv
        tsv = open(tsv_name, "w")
        tsv.write("\t".join(header) + "\n")  # Write the header line
        records = list(SeqIO.parse(gbk, "genbank"))  # Read a GenBank file from the standard input and convert it into a list of SeqRecord objects
        for r in records:  # Each record (r) is a contig with a unique LOCUS name in the GenBank file.
            contig = r.name  # LOCUS name
            for f in r.features:  # Iterate through every feature of the current contig.
                feature_type = f.type
                if feature_type in features:
                    # Fetch the locus_tag
                    if "locus_tag" in f.qualifiers:
                        locus_tag = f.qualifiers["locus_tag"][0]
                    else:
                        locus_tag = "unnamed"
                    
                    # Determine which DNA strand the current feature is located in
                    if f.location.strand == 1:
                        strand = "+"
                    else:
                        strand = "-"
                    
                    # Determine whether the current gene is pesudo
                    if "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
                        is_pesudo = "Y"  # Yes
                    else:
                        is_pesudo = "N"  # No
                    
                    # Determine the product name    
                    if "product" in f.qualifiers:
                        product = f.qualifiers["product"][0]
                    else:
                        product = "unknown"
                    
                    if "gene" in f.qualifiers:
                        gene = f.qualifiers["gene"][0]
                    else:
                        gene = "unknown"
                            
                    # Construct the line to be written into the output file                 
                    line = [contig, locus_tag, f.type, str(f.location.start + 1), str(f.location.end), strand, is_pesudo, product, gene]
                    if args.nucl_seq:
                        line += [str(f.extract(r.seq))]
                    if feature_type == "CDS" and args.prot_seq:
                        if "translation" in f.qualifiers.keys():
                            line += [f.qualifiers["translation"][0]]
                        else:  # Pseudo genes may not have any translation.
                            line += ["unknown"]
                    
                    tsv.write("\t".join(line) + "\n")
        tsv.close()
        
    return


def get_input_filenames(gbks):
    gbk_list = list(gbks)
    if len(gbk_list) == 1 and gbk_list[0].startswith("*"):  # *.gbk
        gbk_list = glob.glob(os.path.join(".", gbk_list[0]))  # Get names of all GenBank files under the current working directory
        
    return(gbk_list)


if __name__ == "__main__":
	main()
