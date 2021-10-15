#!/usr/bin/env python
"""
Creating symbolic links according to a table of three columns: sample names, R1 read files and R2 read files, separated by tab characters.
Warning: existing symbolic links and files will be replaced by new links or files if --update/-u is flagged.
Dependencies: Python 3, bash environment

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 6 Oct 2021; the latest update: 15 Oct 2021
"""

import os
import sys
import subprocess
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "Create symbolic links for paired-end read sets")
    parser.add_argument("--tsv", "-t", dest = "tsv", type = str, required = True, help = "A tab delimited, headerless file of five columns: sample name, input read directory, R1 read file, R2 read file, action (Copy or Link)")
    parser.add_argument("--outdir", "-o", dest = "outdir", type = str, required = "True", help = "Path to the output directory (without the forward slash at the end) in which symbolic links will be created or files will be copied into")
    parser.add_argument("--update", "-u", dest = "update", action = "store_true", help = "Update existing links or override existing files using the new read files. Default: skip existing links or files.")
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
            i, in_dir, r1, r2, action = r.split("\t")
            proceed = True
        except ValueError:
            print(f"Warning: Line {line_num} '{r}' does not contain the five values required. Skip this line.", file = sys.stderr)
            proceed = False
        if proceed:
            r1 = os.path.join(in_dir, r1)
            r2 = os.path.join(in_dir, r2)
            t1 = os.path.join(outdir, i + suffix_r1)
            t2 = os.path.join(outdir, i + suffix_r2)
            if action == "Link":
                readset_num += create_link(t1, r1, update) + create_link(t2, r2, update)
            else:
                readset_num += copy_file(t1, r1, update) + copy_file(t2, r2, update)
    print(f"Symbolic links or files were created for {readset_num} read files.")
    return


def create_link(t, r, u):
    if os.path.exists(r):
        if os.path.exists(t):
            if u:  # Replace existing symbolic links and files with new symbolic links
                if os.path.islink(t):
                    print(f"Updating link {t}")
                    subprocess.run(["unlink", t])
                else:
                    print(f"Warning: {t} is an existing file. It is deleted and a symbolic link is created with {r}.", file = sys.stderr)
                    subprocess.run(["rm", t])
                subprocess.run(["ln", "-s", r, t])
                c = 1
            else:
                print(f"Warning: skip existing link/file {t}.", file = sys.stderr)
                c = 0
        else:
            subprocess.run(["ln", "-s", r, t])
            c = 1
    else:
        print(f"Error: read file {r} does not exist. So link {t} was not created.", file = sys.stderr)
        c = 0
    return c


def copy_file(t, r, u):
    if os.path.exists(r):
        if os.path.exists(t):
            if u:  # Replace existing links and files with copied files
                if os.path.islink(t):
                    print(f"Warning: target {t} is a symbolic link. Remove this link and copy file {r}.", file = sys.stderr)
                    subprocess.run(["unlink", t])
                else:
                    print(f"Warning: overrode existing file {t} with file {r}.", file = sys.stderr)
                    subprocess.run(["rm", t])
                subprocess.run(["cp", r, t])
                c = 1
            else:
                print(f"Warning: skip existing link/file {t}.", file = sys.stderr)
                c = 0
        else:
            subprocess.run(["cp", r, t])
            c = 1
    else:
        print(f"Error: read file {r} does not exist, so it is not copied.", file = sys.stderr)
        c = 0
    return c


if __name__ == "__main__":
    main()