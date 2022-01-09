#!/usr/bin/env python
"""
This script generates summary statistics of nodes for a list of input GFA files
that are produced by the assembler SPAdes. It prints a tab-delimited file to the
STDIN. See gfa-spec.github.io/GFA-spec/GFA1.html for specifications of the
GFA format.

Example command line:
    python gfa_stats.py *.gfa > gfa_stats.tsv
    python gfa_stats.py *.gfa | sed 's/__scaffolds.gfa//g' > gfa_stats.tsv
    python gfa_stats.py *.gfa 1> gfa_stats.tsv 2> gfa_stats.err

This script does not consider the overlap length (e.g., 77M), so every node length
reported by this script includes the overlap length. For singleton nodes (namely,
nodes that do not connect to any other nodes), the overlap length = 0.

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 July 2021; the latest update: 16 July 2021
"""

import os
import sys
import glob
from collections import namedtuple

def main():
    """
    The way each OS deals with the wildcard differs. For example, Win10 passes the string '*.gfa' directly to the script
    as sys.argv[1], whereas Linux replaces this wildcard express with all files globed.
    """
    if len(sys.argv) > 2:
        gfas = sys.argv[1 : ]
    else:
        gfas = glob.glob(sys.argv[1])  # sys.argv[1] = "*.gfa"
    print("Summarise nodes in %i GFA files" % len(gfas), file = sys.stderr)

    # Print summary statistics
    print("\t".join(["Assembly", "Node", "Length", "Depth", "Kmers", "Singleton"]), file = sys.stdout)  # The header line
    for g in gfas:  # Filenames are used for filling the first column.
        if os.path.exists(g):
            summarise_gfa(g)
        else:
            print("Error: skip the inaccessible file " + g, file = sys.stderr)
    return


def summarise_gfa(gfa):
    """ Produces summary statistics for a GFA file """
    nodes = dict()
    linked_nodes = set()
    Node = namedtuple("Node", ["length", "depth", "kmers"])

    # Extracts and transforms the S (segment) and L (link) fields of each GFA file
    g = open(gfa, "r")
    line = g.readline().strip()
    while line:
        if line.startswith("S"):  # A segment line
            _, node_name, seq, depth, kmers = line.split("\t")
            nodes[node_name] = Node(length = str(len(seq)), depth = depth[5 : ], kmers = kmers[5 : ])  # Drop "DP:f:" and "KC:i:"
        elif line.startswith("L"):  # A link line, which always comes after the "S" lines
            _, f, _, t, _, _ = line.split("\t")
            if f != t:
                linked_nodes = linked_nodes.union({f, t})
            else:
                linked_nodes.add(f)
        else:
            break  # The "P" lines make up the last section in the GFA file.
        line = g.readline().strip()
    g.close()

    # Mark non-singleton nodes
    for node_name, node_stats in nodes.items():
        is_singleton = "0" if node_name in linked_nodes else "1"  # "0": no; "1": yes.
        print("\t".join([gfa, node_name, node_stats.length, node_stats.depth, node_stats.kmers, is_singleton]), file= sys.stdout)
    return


if __name__ == "__main__":
    main()
