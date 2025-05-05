#!/usr/bin/env python3

# This script takes a fasta file and counts the lengh of each sequence and returns a csv table of gene name, length
# took it from: https://www.biostars.org/p/118954/
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

FastaFile = open(sys.argv[1], 'r')
LenFile = open(sys.argv[2], 'w')

for name, seq in SimpleFastaParser(FastaFile):
    seqLen = len(seq)
    LenFile.write(name + ',' + str(seqLen) + '\n')

FastaFile.close()
LenFile.close()
