#!/usr/bin/env python3

# INFO: 
# Script to process the mmseqs m8 outpit and calculate the GRR (gene relatedness)

from Bio import SeqIO   
import argparse
import sys 
import pandas as pd
import os

def get_args():
    parser = argparse.ArgumentParser(description='Process mmseqs output and calculate GRR')

    # Process arguments
    parser.add_argument('-f', '--filter_m8',
                        action='store_true',
                        help='Run the filtering function')

    parser.add_argument('-w', '--wGRR',
                        action='store_true',
                        help='Run the wGRR calculation')

    # File arguments
    parser.add_argument('-m',
                        '--m8_file',
                        metavar='FILE',
                        required=True,
                        #type=argparse.FileType('rt'),
                        type=str,
                        help='Table resulting from mmseqs2 bidirectional (.m8)')

    parser.add_argument('-a',
                        '--proteins_file_a',
                        metavar='FILE',
                        required=False,
                        type=argparse.FileType('rt'),
                        help='Fasta file of proteins in group a')
    
    parser.add_argument('-b',
                        '--proteins_file_b',
                        metavar='FILE',
                        required=False,
                        type=argparse.FileType('rt'),
                        help='Fasta file of proteins in group b')

    parser.add_argument('-o',
                        '--out_file',
                        help='Output file name',
                        required=False,
                        metavar='FILE',
                        default='split',
                        type=str)

    return parser.parse_args()

# # --------------------------------------------------
def main():

    # initialize command line arguments
    args = get_args()
    
    # Filter the m8 file
    if args.filter_m8:
                    filtered_m8 = process_protein_matches(args.m8_file)
                    print(f"Kept {len(filtered_m8)} best matches after filtering")
                    print(filtered_m8)
                    filtered_m8.to_csv(args.out_file, sep='\t', index=False)
    else:
                    print("Filtering  skipped. Use --filter_m8 flag to run the filtering.")
                    



# # --------------------------------------------------
def process_protein_matches(m8_file):
    """
    Process a protein match file to:
    1. Remove self-matches (query == target)
    2. Keep only the best match for each query (lowest e-value, then highest fident)

    Args:
        file_path: Path to the input file

    Returns:
        DataFrame with filtered results

    """

    # Define column names
    columns = ["query", "target", "qcov", "tcov", "fident", "evalue", "bits"]


    # Read the entire file
    print(f"Reading file: {m8_file}")
    df = pd.read_csv(m8_file, sep='\t', names=columns)
    print(f"Total rows: {len(df)}")

    # Remove self-matches
    filtered_df = df[df['query'] != df['target']]
    print(f"Rows after removing self-matches: {len(filtered_df)}")
    
    # Sort by query, evalue (ascending), and fident (descending)
    sorted_df = filtered_df.sort_values(by=['query', 'evalue', 'fident'], ascending=[True, True, False])
    
    # Keep only the first occurrence of each query (which will be the best match)
    best_matches = sorted_df.drop_duplicates(subset=['query'], keep='first')

    return best_matches

# # --------------------------------------------------   
if __name__ == "__main__":
    main()
