#!/usr/bin/env python

"""
Remove descriptions from sequence headers in a FASTA file.

Example command:

cat seq.fna | python rmSeqDescr.py > seq_out.fna

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 15 June 2020; the latest modification: 15 June 2020
"""

import sys
import re

def main():
    for line in sys.stdin:
        line = line.rstrip()
        if line.startswith(">"):
            line, _ = line.split(" ", 1)  # Split the header by the first space character.
        print(line)

    return


if __name__ == "__main__":
    main()
