# script to parse the results of prodigal to the gene-to-genome file for vContact2
# Mery Touceda Suarez
# May 1st 2023
# desired outcome: a gene to genome file with protein_id, contig_id, keywords (sample)


# load libraries 
library(tidyverse)
library(seqinr) 
library(stringr)


# make this script executable on command line
args<-commandArgs(TRUE)

## Input files and directories
prodigal_output_file=args[1]
name=args[2]

# read files to objects
prodigal_output <- read.fasta(prodigal_output_file) 

# extract the names of the genes in the fasta of prodigal output
protein_id <- names(prodigal_output)
# extract a list of the contigs from the genes
samples <- lapply(protein_id, function(x) str_extract(x,"^[^.]+"))
contigs <- lapply(protein_id, function(x) str_extract(x,"[^_]*_[^_]*_[^_]*"))

# bind that
genes_df <- as.data.frame(protein_id)
genes_df$contig_id <- contigs
genes_df$keywords <- samples
genes_df <- apply(genes_df,2,as.character)

# write result to csv
outfile <- paste(name, "gene-to-genome.csv", sep="_")
write.csv(genes_df, outfile, quote = FALSE, row.names = FALSE)