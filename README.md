# BINF_toolkit
Small tools for bioinformatics developed by Yu Wan.

===

List of scripts:
    gc.py: calculates the length, GC content, and entropy for each record in a multi-fasta file.
    pull_nucl_region.py: extracts a region of nucleotides by positions from a fasta file.
    get_gene_seq.py: extracts gene sequences from a GenBank file, in accordance with a list of (locus_tag, feature type) tuples.
    gbk2tbl.py: converts a GenBank file into a Sequin feature table that can be processed by tbl2asn used for NCBI WGS submission.