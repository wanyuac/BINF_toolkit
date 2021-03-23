#!/bin/bash

# Copyright (C) 2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Release: 23 Mar 2021

display_usage(){
    echo "
    Assuming SPAdes have produced assembly files for each genome under a subdirectory named by the genome.
    For example, for genomes g1, g2, and g3, subdirectories g1/, g2/, and g3/ have been created under a
    parental directory assemblies/.
    
    Command line: sh saveSPAdesAssemblies.sh [parental directory]
    Example command line: sh saveSPAdesAssemblies.sh \$PWD

    Output: essential assembly files will be copied from subdirectories to the parental directory and be renamed
    by subdirectory names (genome names). Then users may delete subdirectories to save space.
    "
}

d=$1  # Parental (output) directory
cd $d
for dsub in `ls -1 -d */`; do
    g=`basename $dsub`  # Remove the end '/' character and use the subdirectory name as the genome name

    # Assembly graphs
    cp $g/assembly_graph.fastg ./${g}.fastg
    cp $g/assembly_graph_after_simplification.gfa ./${g}__simplified.gfa
    cp $g/assembly_graph_with_scaffolds.gfa ./${g}__scaffolds.gfa  # May be the same as assembly_graph_after_simplification.gfa.

    # Contigs
    cp $g/contigs.fasta ./${g}__contigs.fna  # It has less nodes than does assembly_graph.fastg.
    cp $g/contigs.paths ./${g}__contigs.paths

    # Scaffolds
    cp $g/scaffolds.fasta ./${g}__scaffolds.fna  # It has less contigs (scaffolded) than does contigs.fasta.
    cp $g/scaffolds.paths ./${g}__scaffolds.paths

    # Supplementary information
    cp $g/spades.log ./${g}.log
done
