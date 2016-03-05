#!/usr/bin/python
# Author: Yu Wan, University of Melbourne, 2014-3-21~4-1, 12 May 2015
# Contact: wanyuac@126.com
# This program calculates the length, GC content, and entropy for each record in a multi-fasta file.
# GitHub: https://github.com/wanyuac/BINF_toolkit
# Input: a fasta file which contains multiple sequences from the standard input
# Output: for each sequence, print: 1) header 2) total sequence length 3) percentage of G+C 4) entropy of the sequence
# Command line: python gc.py < filename.fasta
# Treatment of the extended alphabet:
#    1) consider all of 15 characters
#    2) construct a weighted-count table using dictionary
#    3) for each character in the table, take the probability of being A, G, C or T as effective counts
#    4) counts for A, G, C and T is computed by adding up the vectors for every character read from the sequence.
# Licence: GNU GENERAL PUBLIC LICENSE 2.0

import sys
import math

alphabet = {  # weighted-count table
#                 A  G  C  T
            'A': [1, 0, 0, 0],
            'G': [0, 1, 0, 0],
            'C': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'S': [0, 0.5, 0.5, 0],
            'W': [0.5, 0, 0, 0.5],
            'R': [0.5, 0.5, 0, 0],
            'Y': [0, 0, 0.5, 0.5],
            'M': [0.5, 0, 0.5, 0],
            'K': [0, 0.5, 0, 0.5],
            'V': [0.33, 0.33, 0.33, 0],
            'H': [0.33, 0, 0.33, 0.33],
            'D': [0.33, 0.33, 0, 0.33],
            'B': [0, 0.33, 0.33, 0.33],
            'N': [0.25, 0.25, 0.25, 0.25]
            }

def read_fasta (fasta):
    header = []  # starts from 0
    s = 'start'
    seq = []
    for line in fasta:
        line = line.rstrip('\n')  # remove '\n' at the end
        if line.startswith('>'):  # find a new sequence
            header.append(line)  # Do not use header = ... here because .append() returns None
            seq.append(s.upper())  # Append the last row of last sequence, note that the first element is ''.
            s = ''  # reset s
        else:
            s = s + line.upper()
    seq.append(s)  # append the last string to list seq after the loop
    seq = seq[1:]  # remove the first element
    list = [header, seq]
    return list

def base_count (seq):
    #      A  G  C  T
    num = [0, 0, 0, 0]  # numbers of A, G, C, T
    L = len(seq)
    for i in range(0, L):
        num = [sum(j) for j in zip(num, alphabet[seq[i]])]  # addition of vectors: element by element
    return num

def GC_content (num, seq_len):
    GC = float(num[1] + num[2])  # sum of the numbers of G and C in the list num.
    return GC / seq_len

def entropy (num, seq_len):
#        P(A) P(G) P(C) P(T)
    p = [0, 0, 0, 0]  # initiate a list
    H = 0  # entropy
    y = 0
    for i in range(0, 4):
        p[i] = float(num[i]) / seq_len  # estimates the probability by frequency
        if p[i] == 0:  # cannot take log0
            y = 0  # lim(x^x) = 0 when x -> 0
        else:
            y = p[i] * math.log(p[i], 2)
        H = H - y
    return H

#/////////////// Main program /////////////////////
# read from the standard input
content = sys.stdin.readlines()  # including '\n'
fasta = read_fasta(content)  # fasta[0]: headers, fasta[1]: sequences
header = fasta[0]
seq = fasta[1]
n = range(0, len(header))  # number of sequences
for i in n:
    s = seq[i]
    L = len(s)
    num = base_count(s)
    print header[i]  # output 1: headers
    print L  # output 2: the length of the sequence
    print '%04.2f'%(GC_content(num, L) * 100)  # output 3: the G+C content of this sequence
    print '%02.1f'%entropy(num, L)  # output 4: the entropy