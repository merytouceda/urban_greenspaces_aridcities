#!/usr/bin/env python3

# INFO: 
# Script to retrieve the fasta sequences belonging to contigs that I have annotated as plasmids

from Bio import SeqIO   
import argparse
import sys 
import pandas as pd
import os

def get_args():
    parser = argparse.ArgumentParser(description='Retrieve fasta sequence of contigs that are plasmids')

    parser.add_argument('-p',
                        '--plasmids_file',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help='List of contigs that are plasmids (.csv)')


    parser.add_argument('-c',
                        '--contigs_file',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help='Fasta file of contigs')

    parser.add_argument('-o',
                        '--out_dir',
                        help='Output directory',
                        metavar='DIR',
                        default='split',
                        type=str)

    return parser.parse_args()

# # --------------------------------------------------
def main():

    # initialize command line arguments
    args = get_args()


    # 1. Create dictionaries

    # Contigs dictionary (from fasta, see function below)
    contigs_dict = fasta_to_dict(args.contigs_file)

    # Plasmids dictionary (from .csv)
    plasmidsFile=pd.read_table(args.plasmids_file, header = 0, sep='\t')
    plasmids_dict = dict(zip(plasmidsFile.Contig,plasmidsFile.x)) # (i don't really care about x, it is just the row number)

    # Create output parameters
    basename = os.path.basename(args.contigs_file.name)
    root, ext = os.path.splitext(basename)

    # Find plasmids and mges contigs in contig fasta file and write the ones that are found to output fasta file
    with open(os.path.join(args.out_dir, root + 'mobile_contigs' + ext),
                "wt", encoding='utf-8') as outFile:
        for p_key in plasmids_dict:  # iterate over the list of contig headers
            for c_key in contigs_dict:
                if p_key == c_key: 
                    outFile.write(str(">" + c_key + '\n'))
                    outFile.write(str(contigs_dict[c_key] + '\n'))
        
        


# # --------------------------------------------------
def fasta_to_dict(fasta_file):
    """
    Reads a FASTA file and returns a dictionary where keys are sequence IDs 
    and values are the corresponding sequences.
    """

    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        modified_header = record.id.replace('.contigs.fa', '') # fix the contig name in the fasta file
        modified_header = modified_header.replace('./', '')
        modified_header = modified_header.replace('flag(.*)', '')     
        seq_dict[modified_header] = record.seq

    return seq_dict


# # --------------------------------------------------   
if __name__ == "__main__":main()