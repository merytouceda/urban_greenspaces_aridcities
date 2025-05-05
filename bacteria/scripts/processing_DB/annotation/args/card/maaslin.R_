# This script runs maaslin2 for urban dataset

library(Maaslin2)

args<-commandArgs(TRUE)

## Input files and directories
metadata=read.csv(args[1], row.names = 1) # full path to metadata
count_table=read.csv(args[2], row.names = 1) # full path to count table

fit_data = Maaslin2(
  input_data = count_table, 
  input_metadata = metadata, 
  output = "carbon_maaslin", 
  fixed_effects = "LaNdUse4.1", 
  random_effects = "Paired_ID_v4")