"""
This script screens genes against a reference database using megaBLAST for every input FASTA file.
Specifically, it takes as input a list of FASTA files and searches every DNA sequence against the reference database. Obviously,
it performs a targeted analysis for input sets of DNA sequences.
For example, you may want to find out all resistance genes in every bacterial genome, for which you may create FASTA files of
coding sequences for every genome and use this script to profile this kind of genes.

Number of options: 6 (2 compulsory and 4 optional)
Usage:
    python screen_genes_blast.py --in *.fna --db [reference database] --strains [a comma-delimited string of strain names]
    --genomes [a comma-delimited string of genome names] --opt [options and arguments for BLAST]
    --outfmt [output format code] > [output file name]
    
Prerequisite: A BLAST nucleotide database should be made before using this script.
    makeblastdb -in your.fasta -dbtype nucl -out db_name -logfile your.log

Options "--strains" and "--genomes" are optional.

A spreadsheet can be created beforehand to ensure the strain name and the genome name to match each FASTA file:
    strain  genome  fasta_file
    AH0650_Sm1  chr chr.fna
    AH0650_Sm1  plasmid plasmid.fna

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Development history: 3 July 2016
Python version: 2.7.10
License: GNU GPL 2.1
"""

from argparse import ArgumentParser
import sys, os, subprocess

def parse_arguments():
    # read arguments of options
    parser = ArgumentParser(description="Fix problems in SRST2's ARG-Annot database")
    parser.add_argument("--in", "-i", dest = "input", nargs='+', type = str, required = True, default = "", help = "A list of input FASTA files")
    parser.add_argument("--db", "-d", dest = "db", type = str, required = True, default = "", help="A reference nucleotide database for BLAST")
    parser.add_argument("--strains", "-s", dest = "strains", type = str, required = False, default = "", help = "(optional) Comma-delimited names of bacterial strains")
    parser.add_argument("--genomes", "-g", dest = "genomes", type = str, required = False, default = "", help = "(optional) Comma-delimited genome names")
    parser.add_argument("--opt", "-o", dest = "opt", type = str, required = False, default = "-evalue 0.001 -max_target_seqs 2 -perc_identity 98",\
                        help = "Options and argument passed to BLAST")
    parser.add_argument("--outfmt", "-f", dest = "outfmt", type = str, required = False,\
                        default = "6 qseqid sseqid qstart qend sstart send qlen slen length bitscore pident qcovs gaps evalue",\
                        help = "The configuration of the 'outfmt' option for BLAST")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    n_fasta = len(args.input)
    
    # parse strain information
    if args.strains != "":
        strains = args.strains.split(",")
        n_str = len(strains)
    else:
        strains = None
        n_str = 0
        
    # parse genome information
    if args.genomes != "":
        genomes = args.genomes.split(",")
        n_gen = len(genomes)
    else:
        genomes = None
        n_gen = 0
    
    # check whether strains, genomes and files match
    if  n_str != n_fasta:
        sys.exit("Error: strain number is not equal to the number of FASTA files.")
        
    if n_gen != n_fasta:
        sys.exit("Error: genome number is not equal to the number of FASTA files.")
    
    # get column names of the output file
    colnames = args.outfmt.split(" ")[1 : ]  # remove the first element -- the format id
    
    # print the header line to the stdout
    if n_gen > 0:
        colnames = ["genome"] + colnames
    if n_str > 0:
        colnames = ["strain"] + colnames
    print "\t".join(colnames)
    
    # search every set of query sequences against the reference database
    i = 0  # the counter of FASTA files
    for fasta in args.input:
        cmd = ["blastn", "-task", "megablast", "-db", args.db, "-query", fasta] + \
               args.opt.split(" ") + ["-outfmt", args.outfmt]  # Each pair of the option and its argument must be separated as elements of a list.
        proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out = proc.communicate()  # obtain the output of BLAST from the standard output
        hits = out[0].splitlines()  # stderr: out[1]
        
        # print all lines in the current output
        for line in hits:
            if n_gen > 0:
                line = genomes[i] + "\t" + line  # add the genome name to each line
            if n_str > 0:
                line = strains[i] + "\t" + line  # add the strain name to each line
            print line
        i += 1
    
if __name__ == "__main__":
    main()
    