#!/usr/bin/env python3

# INFO: 
# Script to process the mmseqs m8 output and calculate the GRR (gene relatedness)

from Bio import SeqIO   
import argparse
import sys 
import pandas as pd
import os

def get_args():
    parser = argparse.ArgumentParser(description='Process mmseqs output and calculate GRR')

    # Process arguments
    parser.add_argument('-fm', '--filter_mixed',
                        action='store_true',
                        help='Run the extraction of only mixed hits (different group in query and target)')

    parser.add_argument('-fb', '--filter_bbh',
                        action='store_true',
                        help='Run extraction of only best bi-directional hits')

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
                        help='Filtered table resulting from mmseqs2 bidirectional (.m8)')

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

    # Open the output file in **append mode** so we don’t overwrite existing content
    with open(args.out_file, "a") as f:
        
        # Run the mixed filter
        if args.filter_mixed: 
            mixed_filtered_m8_file = filter_mixed_types_only(args.m8_file)
            print(f"Kept {len(mixed_filtered_m8_file)} mixed matches after filtering")
            print(mixed_filtered_m8_file)
        else: 
            mixed_filtered_m8_file = pd.read_csv(args.m8_file, sep='\s+') # just read input file into a DataFrame if not filtering
            print("Skipping the extraction of different group matches (different group query and target). Use --filter_mixed flag to run the filtering.")

        # Run the BBH filter
        if args.filter_bbh:
            bbh_m8_file = extract_bbh(mixed_filtered_m8_file)
            print(f"Kept {len(bbh_m8_file)} BBHs  after filtering")
            print(bbh_m8_file)
        else:
            bbh_m8_file = mixed_filtered_m8_file
            print("Skipping the extraction of best bi-directional hits (BBH). Use --filter_mixed flag to run the filtering.")

        # Calculate the wGRR
        if args.wGRR:
            wGRR = calculate_wGRR(bbh_m8_file, args.proteins_file_a, args.proteins_file_b)
            print(wGRR)
            
            # Save the m8_file name and wGRR result to the output file
            f.write(f"{args.m8_file}\n")  # Write input file name
            f.write(f"{wGRR}\n")  # Write wGRR result

        # Calculate Sorensen

        # Calculate Jaccard
    
                    
# # -------------------------------------------------- 

def filter_mixed_types_only(filtered_m8_file):
    """
    Filter the DataFrame to keep only rows where:
    - One field is a contig and the other is a plasmid
    - Removes rows where both are contigs or both are plasmids
    
    Args:
        filtered_m8_file: filtered m8 table with "query" and "target" columns                                                                                                                                                                                              
        NOTE: in case of having done the filtering on chunks, make sure you concatenate them properly with:                                                                                                                                                                  
        awk 'FNR==1 && NR!=1{next}{print}' _filtered.txt > combined_file.txt

    Returns:
        Filtered pandas DataFrame
    """
    # Check if input is a DataFrame or a file path
    if isinstance(filtered_m8_file, pd.DataFrame):
        df = filtered_m8_file
    else:
        df = pd.read_csv(filtered_m8_file, sep='\s+')

    # Create helper functions to check if a string contains "contig" or "plasmid"
    def is_contig(s):
        return "CONTIG" in s.upper()
    
    def is_plasmid(s):
        return "PLASMID" in s.upper()
    
    # Create masks for different combinations
    both_contigs = df['query'].apply(is_contig) & df['target'].apply(is_contig)
    both_plasmids = df['query'].apply(is_plasmid) & df['target'].apply(is_plasmid)
    
    # Keep only rows where not (both are contigs or both are plasmids)
    mixed_types = ~(both_contigs | both_plasmids)
    
    return df[mixed_types].reset_index(drop=True)

# # --------------------------------------------------                                                                                                                                                                                                                        
def extract_bbh(mixed_filtered_m8_file):
    """                                                                                                                                                                                                                                                                       
    Extract bi-directional best hits (BBH) from a filtered m8 table

    Args: 
       mixed_filtered_m8_file: filtered m8 file of only mixed (query and target from different groups/organisms)


    Returns: 
        DataFrame of best bidirectional hits (BBH)

    """

    # Check if input is a DataFrame or a file path
    if isinstance(mixed_filtered_m8_file, pd.DataFrame):
        df = mixed_filtered_m8_file
    else:
        df = pd.read_csv(mixed_filtered_m8_file, sep='\s+')
    
    # Create a set of all query-target pairs
    pairs = set(zip(df['query'], df['target']))
    
    # Find reciprocal matches
    reciprocal_matches = []
    for query, target in pairs:
        # Check if the reverse pair exists
        if (target, query) in pairs:
            # Find the row indices for both directions
            forward_match = df[(df['query'] == query) & (df['target'] == target)]
            reverse_match = df[(df['query'] == target) & (df['target'] == query)]
            
            # Add one of them to results if both exist
            if not forward_match.empty and not reverse_match.empty:
                reciprocal_matches.append(forward_match.iloc[0])
                #reciprocal_matches.append(reverse_match.iloc[0])
    
    # Create a DataFrame from the matches
    if reciprocal_matches:
        result_df = pd.DataFrame(reciprocal_matches)
        return result_df
    else:
        return pd.DataFrame(columns=df.columns)    


# # --------------------------------------------------
def calculate_wGRR(bbh_m8_file, a, b):
    """
    Process a protein match file to calculate the wGRR (weighted gene repertoire relatedness)

    wGRR(A,B) = sum(pident(A,B))/min(len(A),len(B))
    
    where: 
    A,B = proteins in group/organism A and B respectively
    pident(A,B) = percent identity between A and B (column in m8 file);
    NOTE: we will only use the A,B pairs that are the best bidirectional hit (the best match of eachother)
    len(A), len(B) = the number of proteins in group A and B
 
    Args:
        bbh_m8_file: Filtered m8 table of best bidirectional hits
        NOTE: this file should be a concatenation of the filtered chunks of original m8 file
        that has gone through the bbh_extraction function
        a: fasta file for group a (ex. plasmid proteins fasta)
        b: fasta file for group b (ex. bacterial proteins fasta)
        NOTE: these should be the fasta files used to create the m8 table

    Returns:
        a number = wGRR

    """
    # open m8 file
    if isinstance(bbh_m8_file, pd.DataFrame):
        m8_df = bbh_m8_file
    else:
        m8_df = pd.read_csv(bbh_m8_file, sep='\s+')

    # open a fasta file and calculate length
    num_prots_a = sum(1 for _ in SeqIO.parse(a, "fasta"))    
    # open b file and calculate length
    num_prots_b = sum(1 for _ in SeqIO.parse(b, "fasta"))

    # calculate wGRR
    fident_sum = m8_df["fident"].sum()
    wGRR = fident_sum/min(num_prots_a, num_prots_b)

    return(wGRR)

# # --------------------------------------------------   
if __name__ == "__main__":
    main()
