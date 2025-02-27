#!/usr/bin/env python3

from Bio import SeqIO   
import argparse
import sys                                                                                                                                                                                                                   


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Extract head and tail of fasta sequences from multifasta file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_file',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help='Input multi-fasta file')

    parser.add_argument('-n',
                        '--number',
                        help='Proportion of bp to extract',
                        type=int,
                        default=0.25)

    parser.add_argument('-oh',
                        '--outfile_head',
                        help='Output file for the head section of sequences (first n basepairs)',
                        metavar='FILE',
                        type=argparse.FileType('wt'))

    parser.add_argument('-ot',
                        '--outfile_tail',
                        help='Output file for the tail section of sequences (last n basepairs)',
                        metavar='FILE',
                        type=argparse.FileType('wt'))

    return parser.parse_args()

# --------------------------------------------------
def main():

    """From each fasta sequence in file retrieve the first and last n basepairs
    save each one in a file"""

    args = get_args()

    for seq_record in SeqIO.parse(args.input_file, "fasta"):
        # get the sequence from multifasta entry (fasta sequence)
        sequence = str(seq_record.seq)

        # stablish the number of basepairs to select as a proportion of the sequence's length
        new_n = int(round(len(sequence)*args.number, 0))

        # save the first and last n basepairs in an object
        head_sequence = sequence[:new_n]
        tail_sequence = sequence[-new_n:]

        # save sequence and record multifasta for each
        head_fasta_record = f'>{seq_record.id}\n{head_sequence}\n'
        tail_fasta_record = f'>{seq_record.id}\n{tail_sequence}\n'

        # write formed fasta record to respective output file
        args.outfile_head.write(head_fasta_record)
        args.outfile_tail.write(tail_fasta_record)

    print(f'All done!')


                    
                                                                                                                                                                                                                            
# --------------------------------------------------
if __name__ == '__main__':
    main()