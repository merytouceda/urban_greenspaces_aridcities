# This script takes the input directory with all the balst.out results and merges them into one file and creates a count table of carbon genes in the samples
# This script is specifically made for the CARBON data in the local urban (LANDUSE, LU) dataset

# load libraries
library(tidyverse)

args<-commandArgs(TRUE)

## Input files and directories
input_dir=args[1] # directory with all blast output files
count_table=args[2] # full path to output count table


setwd(input_dir)
# Read all the files and save them in one big file
# add column names based on blast fmt6
colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
              "qend", "sstart", "send", "evalue", "bitscore")



# do it with all the files: 

# initialize a data frame to add all others and make a huge one
colnames2 <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
               "qend", "sstart", "send", "evalue", "bitscore", "sample_id")
d0 <- data.frame(matrix(NA, nrow = 1, ncol = 13))
colnames(d0) <- colnames2

# list of files: 
files <- list.files(pattern="*_blastx_strict.out", full.names=F, recursive=FALSE)

# make a full table of all the samples
for (file in files){
  name <- basename(file)
  d <- tryCatch(read.table(file, header = F, sep = "\t"), error=function(e) NULL)
  colnames(d) <- colnames
  # extract sample name from name of file
  samplename <- str_extract(name, "[^_]*")
  
  d <- d %>%
    group_by(qseqid) %>% # group by same read name
    arrange(evalue,.by_group = T) %>% # sort by evalue within group
    slice_head(n = 1) %>% # select the first hit (best one) for each group
    mutate(sample_id=samplename)
  
  # append to form a big dataset
  d0 <-rbind(d0, d)
  
}

print(head(d0))

# get the count table
cazy <- d0 %>%
  filter(!is.na(sseqid)) %>%
  group_by(sseqid,sample_id) %>%
  tally() %>%
  pivot_wider(names_from = sample_id, values_from = n, values_fill = 0)%>%
  column_to_rownames(var = "sseqid") 

cazy <- cazy %>%
  select(order(colnames(cazy)))

#change - for . in colnames (so they match metadata)
names(cazy) <- gsub(x = names(cazy), pattern = "-", replacement = ".")

print(head(cazy))


# save output
write.csv(cazy, count_table)