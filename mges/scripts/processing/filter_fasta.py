#!/usr/bin/env python3

import sys
from Bio import SeqIO

def filter_fasta_files(plasmid_file, all_genes_file, output_file):
    # Read plasmid sequences into a set
    plasmid_sequences = set()
    for record in SeqIO.parse(plasmid_file, "fasta"):
        plasmid_sequences.add(str(record.seq))
    
    # Filter all_genes file and write non-matching sequences to output
    with open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(all_genes_file, "fasta"):
            if str(record.seq) not in plasmid_sequences:
                SeqIO.write(record, output_handle, "fasta")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <plasmid_fasta> <all_genes_fasta> <output_fasta>")
        sys.exit(1)
    
    plasmid_file = sys.argv[1]
    all_genes_file = sys.argv[2]
    output_file = sys.argv[3]
    
    filter_fasta_files(plasmid_file, all_genes_file, output_file)

if __name__ == "__main__":
    main()
