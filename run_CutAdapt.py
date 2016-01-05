'''
This script runs CutAdapt for a list of paired-end readsets. Author: Yu Wan (wanyuac@gmail.com) License: GNU GPL 2.1 Edition 
Edition history: 15 Dec 2015, 1-2 Jan 2016
Licence: GNU GPL 2.0
'''

import os, re
from argparse import ArgumentParser

MEMORY = "2048" # 2 GB for each job 
WALL_TIME = "1-0:0:0"

def parse_args():
    parser = ArgumentParser(description= "Run CutAdapt to remove adapter sequences from reads.")
    parser.add_argument("--reads", type = str, required = True, help = "A list of paired-end readsets")
    parser.add_argument("--f_adapter", type = str, required = True, help = "Adapter sequences of forward reads")
    parser.add_argument("--r_adapter", type = str, required = True, help = "Adapter sequences of reverse reads")
    parser.add_argument("--pattern", type = str, required = False, default = "\d\d\d\d_\d#\d*", help = "A regular expression for pulling out sample names")
    parser.add_argument("--side", type = str, required = False, default = "3'", help = "3'-end adapters or 5'-end adapters?")
    parser.add_argument("--len", type = str, required = False, default = "108", help = "Minimun read length")
    parser.add_argument("--overlap", type = str, required = False, default = "33", help = "Minimun overlap length between an adapter and a read")
    parser.add_argument("--discrep", type = str, required = False, default = "0.03", help = "The discrepancy rate between the reference adapter sequences and subjects")
    parser.add_argument("--outdir", type = str, required = True, default = ".", help = "The directory for outputs")
    return parser.parse_args()

def submit_jobs(readsets, adapt_f, adapt_r, side, outdir, min_read_len, discrep, overlap):
    for sample, reads in readsets.iteritems():
        cmd = '#!/bin/bash'
        cmd += '\n#SBATCH -p main'
        cmd += '\n#SBATCH --job-name=CutAdapt'
        cmd += '\n#SBATCH --ntasks=1'
        cmd += '\n#SBATCH --mem-per-cpu=' + MEMORY
        cmd += '\n#SBATCH --time=' + WALL_TIME
        cmd += '\ncd ' + outdir + '\n'
        if side == "3'":
            cutadapt_cmd = 'cutadapt -a file:' + adapt_f + ' -A file:' + adapt_r + ' --minimum-length ' + min_read_len + ' -e ' + discrep + \
                           ' --overlap ' + overlap + ' -o ' + sample + '_1.fastq.gz' + ' -p ' + sample + '_2.fastq.gz' + ' ' + reads[0] + ' ' + reads[1]
        else:
            cutadapt_cmd = 'cutadapt -g file:' + adapt_f + ' -G file:' + adapt_r + ' --minimum-length ' + min_read_len + ' -e ' + discrep + \
                           ' --overlap ' + overlap + ' -o ' + sample + '_1.fastq.gz' + ' -p ' + sample + '_2.fastq.gz' + ' ' + reads[0] + ' ' + reads[1]
        cmd += cutadapt_cmd + ' > ' + sample + '.log'
        #print cmd
        print cutadapt_cmd
        os.system("echo '" + cmd + "' | sbatch")
              
def load_reads(f, pattern):
    readsets = {}
    samples = []
    with open(f, "rU") as inputs:
        lines = inputs.read().splitlines()
        
    # initialises the dictionary
    for line in lines:
        samples.append(re.findall(pattern, line)[0])
    samples = sorted(list(set(samples)))  # remove redundancy and sort the list
    for sample in samples:
        readsets[sample] = ["", ""]
    
    # matches paired-end read sets to samples
    for line in lines:
        sample = re.findall(pattern, line)[0]
        orentation = int(re.findall("_\d\.", line)[0][1])  # e.g. "_1." => 1
        readsets[sample][orentation - 1] = line
        
    return(readsets)
    
def main():
    args = parse_args()
    readsets = load_reads(args.reads, args.pattern) # generates a dictionary {sample:[read1, read2],...}
    submit_jobs(readsets, args.f_adapter, args.r_adapter, args.side, args.outdir, args.len, args.discrep, args.overlap) # submit a job for each pair of readsets

if __name__ == '__main__':
    main()
