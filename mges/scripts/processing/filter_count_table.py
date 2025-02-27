#!/usr/bin/env python3
 
import argparse
import sys
import pandas as pd                                                                                                                                                                                                 


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Extract head and tail of fasta sequences from multifasta file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_table',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help='Input count table')

    parser.add_argument('-f',
                        '--table_seqs_keep',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help='Table of final plasmid sequences to keep')

    parser.add_argument('-oh',
                        '--out_file',
                        help='Output, filtered count table',
                        metavar='FILE',
                        type=argparse.FileType('wt'))

    return parser.parse_args()

# --------------------------------------------------
def main():

    """ Filter the input count table to only keep those sequences in the final plasmid table """

    args = get_args()

   # Read count table into dictionary
    lines = args.input_table.readlines()

    # Assuming the table is space-separated
    header = lines[0].strip().split()
    data_dict = {col: [] for col in header}

    for line in lines[1:]:
        row = line.strip().split()
        for col, value in zip(header, row):
            data_dict[col].append(value)

    print(data_dict)


    # Read plasmid csv table into dictionary
    lines = args.table_seqs_keep.readlines()

    header = lines[0].strip().split(',')
    plasmid_dict = {}

    # Process each row and use the first column as the key
    for line in lines[1:]:
        row = line.strip().split(',')
        plasmid_dict[row[0]] = row[1:] 

    print(plasmid_dict)
 

                                                                                                                                                                                                                                    
# --------------------------------------------------
if __name__ == '__main__':
    main()