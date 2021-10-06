#!/usr/bin/env python
"""
Creating symbolic links according to a table of three columns: sample names, R1 read files and R2 read files, separated by tab characters.
Dependencies: Python 3, bash environment
Copyright (C) 2017-2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
First edition: 6 Oct 2021; the latest update: 6 Oct 2021
"""

import os
import sys
import subprocess
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "Create symbolic links for paired-end read sets")
    parser.add_argument("--tsv", "-t", dest = "tsv", type = str, required = True, help = "A tab delimited, headerless file of three columns: sample name, R1 read file, R2 read file")
    parser.add_argument("--outdir", "-o", dest = "outdir", type = str, required = "True", help = "Path to the output directory (without the forward slash at the end) in which symbolic links will be created")
    parser.add_argument("--update", "-u", dest = "update", action = "store_true", help = "Update existing links using the new read files. Default: skip existing links.")
    parser.add_argument("--R", "-R", dest = "R", action = "store_true", help = "Create symbolic links with suffices _R[1,2].fastq.gz rather than the default _[1,2].fastq.gz")
    return parser.parse_args()


def main():
    params = parse_arguments()

    with open(params.tsv, "r") as f:
        readsets = f.read().splitlines()

    outdir = params.outdir
    update = params.update

    if params.R:
        suffix_r1 = "_R1.fastq.gz"
        suffix_r2 = "_R2.fastq.gz"
    else:
        suffix_r1 = "_1.fastq.gz"
        suffix_r2 = "_2.fastq.gz"
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    line_num = 0
    readset_num = 0
    for r in readsets:
        line_num += 1
        try:
            i, r1, r2 = r.split("\t")
            proceed = True
        except ValueError:
            print(f"Warning: Line {line_num} '{r}' does not contain the three values required. Skip this line.", file = sys.stderr)
            proceed = False
        if proceed:
            t1 = os.path.join(outdir, i + suffix_r1)
            t2 = os.path.join(outdir, i + suffix_r2)
            create_link(t1, r1, update)
            create_link(t2, r2, update)
            readset_num += 1
    print(f"Symbolic links were created for {readset_num} readsets.")
    return


def create_link(t, r, u):
    if os.path.exists(r):
        if os.path.exists(t):
            if u:
                if os.path.islink(t):
                    print(f"Updating link {t}")
                    subprocess.run(["unlink", t])
                    subprocess.run(["ln", "-s", r, t])
                else:
                    print(f"Warning: symbolic link {t} was not created as it is an existing file.", file = sys.stderr)
            else:
                print(f"Skip existing link/file {t}")
        else:
            subprocess.run(["ln", "-s", r, t])
    else:
        print(f"Error: read file {r} does not exist. So link {t} was not created.")
    return


if __name__ == "__main__":
    main()