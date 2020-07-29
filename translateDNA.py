#!/usr/bin/env python3

"""
Translate DNA sequences into protein sequences using a given codon table. The translation terminates at the
first stop codon.

Command:
    translateDNA.py -c [codon table number] -k

Example commands:
cat seq.fna | python translateDNA.py > seq.faa  # Print error messages on screen
cat seq.fna | python translateDNA.py -c 11 -k 1>seq.faa 2>seq.err
cat seq.fna | python translateDNA.py 10 | python rmSeqDescr.py > seq.faa

Default codon table number is 11 (for prokaryotic genes).

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 15 June 2020; the latest modification: 29 July 2020
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser("Convert nucleotide sequences into protein sequences")
    parser.add_argument("-c", type = int, required = False, default = 11, action = "store", help = "Codon table (Default: 11)")
    parser.add_argument("-k", action = "store_true", help = "Use sequence IDs as protein names without capitalise the first letter (Default: off)")

    return parser.parse_args()


def main():
    # Initialisation
    args = parse_arguments()
    codon_tab = args.c
    print("Codon table: %i" % codon_tab, file = sys.stderr)

    # Go through sequences from the stdin
    for rec in SeqIO.parse(sys.stdin, "fasta"):
        # Process the sequence description
        seqid = rec.id
        seqid_len = len(seqid)
        descr = "" if seqid_len == len(rec.description) else rec.description[(seqid_len + 1) : ]  # Drop seqid from rec.description
        seqid_new = seqid if args.k else geneName_to_proteinName(seqid)
        
        # Translate the DNA sequence into a protein sequence
        try:
            # Set cds = True to check error: 'First codon is not a start codon'.
            # Set to_stop = False to print "*" that represents stop codons. Otherwise, the asterisk is not printed.
            rec_prot = rec.translate(table = codon_tab, id = seqid_new, description = descr, to_stop = True, cds = False, stop_symbol = "*")
            trans_succ = True
        except KeyError:
            print("Warning: sequence \"%s\" cannot be translated." % rec.description, file = sys.stderr)
            rec_prot = SeqRecord(Seq(''), id = '', name = '', description = '')
            trans_succ = False

        # Print protein sequences
        if trans_succ:
            if descr == "":
                print(">" + rec_prot.id, file = sys.stdout)
            else:
                print(">%s %s" % (rec_prot.id, rec_prot.description), file = sys.stdout)
            print(rec_prot.seq, file = sys.stdout)

    return


def geneName_to_proteinName(g):
    # Converting a gene name to protein name through capitalising the first letter of the gene name
    p = g[0].upper() + g[1 : ]

    return p


if __name__ == "__main__":
    main()
