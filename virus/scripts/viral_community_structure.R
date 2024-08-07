# Viral community structure visualization and analysis

# load libraries
library(ggplot2)
library(vegan)
library(see)
library(dplyr)
library(prabclus)
library(lme4) #mixed effect models
library(lmerTest)
library(agricolae)
library(car)
library(performance)
library(tidyverse)

# initialize palette
pal1 = c("#FFC20A", "#0C7BDC")
pal2 = c("#009E73","#E66100", "#5D3A9B")


####################################################################################
############## SETUP
####################################################################################
setwd("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/virus/data")
counts <- read.table("viral_species_count_table.txt", header = T, row.names = 1, sep = "\t")
counts <- as.data.frame(t(counts))
metadata<- read.csv("../../metadata.csv", row.names = 1)

# turn the variables of interest into factors 
metadata$urban.natural <- as.factor(metadata$urban.natural)
metadata$landuse2 <- as.factor(metadata$landuse2)

# Fix the counts rownames so they match the metadata rownames 
# check
rownames(counts) == rownames(metadata)
# fix
counts <- counts %>%
  rownames_to_column(var = "sample") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "^X", replacement = "")) %>%
  column_to_rownames(var = "sample")
# check again
rownames(counts) == rownames(metadata)


# remove those contigs with length < 1000bp 
checkv <- read.table("final_checkv.tsv", header = T, sep = "\t")
# find the short contigs
short_contigs <- checkv[checkv$contig_length < 1000, ]
# remove them from checkv
checkv <- checkv[!checkv$contig_length < 1000, ]
# remove them from count table
counts <- counts[,!colnames(counts) %in% short_contigs$contig_id]


# Calculate the number of reads that are viral for publication
metadata$sample == rownames(counts)
mean(rowSums(counts)/metadata$total.reads)
sd(rowSums(counts)/metadata$total.reads)
summary(rowSums(counts)/metadata$total.reads)
# IQR = 0.004242 - 0.002382


####################################################################################
############## GENERAL ECOLOGY
####################################################################################

# ALPHA-DIVERSITY
# summary
summary(rowSums(counts))
hist(rowSums(counts))

# Richness
metadata$viral.species.rich.rar <- specnumber(rrarefy(counts, sample=39718))
# Shannon
metadata$viral.species.shan.rar <- diversity(rrarefy(counts, sample=39718), index = "shannon")

# boxplot
ggplot(metadata, aes(x = urban.natural, y = viral.species.rich.rar))+
  geom_boxplot(alpha=0.6, outlier.shape = NA)+
  geom_jitter(aes(color= vegetation_structure), size = 3)+
  xlab(NULL) + 
  ylab("Number of different viral species")+
  scale_fill_manual(values=pal1)+
  #scale_color_manual(values=pal2)+
  scale_color_manual(values = pal2)+
  #scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))


# pointrange
ggplot(data = metadata, aes(x = urban.natural, y = viral.species.rich.rar, color = urban.natural)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=urban.natural, y=viral.species.rich.rar), position=position_jitterdodge(jitter.width=0.1)) +
  ylab("Number of different viral species")+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  theme_bw()+
  theme(strip.text.x = element_blank(), legend.position = "none")
ggsave("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/virus/figures/viral_rich.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# STATS
Anova(lmer(viral.species.rich.rar ~ urban.natural *landuse2 + (1|site) , data = metadata), type = 3)




# BETA-DIVERSITY

# calculate bray curtis dissimilarity matrix
counts.bray <- vegdist(counts, method="bray")

# Generate ordination with NMDS (non-multidimensional scaling)
counts.nmds <- metaMDS(counts.bray, k=2, try = 100)
metadata$Axis01 = counts.nmds$points[,1]
metadata$Axis02 = counts.nmds$points[,2]
counts.nmds$stress #0.094


# ordination plot
ggplot(metadata, aes(Axis01, Axis02))+
  geom_point(aes(color= vegetation_structure), size=3)+
  stat_ellipse(aes(group = urban.natural), linetype = 2)+
  scale_color_manual(values=pal2)+
  theme_classic()+
  theme(legend.position="none", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/virus/figures/viral_commcomp.pdf", device = "pdf", width = 5, height = 5 , units = "in")

# stats: 
adonis2(counts.bray ~ urban.natural*landuse2, data = metadata, permutations = 999, method = "bray")

# BETA-DISPERSION
mod <- betadisper(counts.bray, metadata$urban.natural)

# mod plot
plot(mod)
# boxplot
boxplot(mod)

# STATS
anova(mod)




########################################
# ABUNDANCE MEASURES
########################################

# RPKM NORMALIZATION
#######       #######      #######     ########     ######   ########    ########
# 1. pivot longer the count table to make every gene in every sample a row
counts_long <- as.data.frame(counts) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = -c("sample"), values_to = "count", names_to= "virus")

# 2. Create table of total reads per sample and load viral length data
total_reads_table <- metadata %>%
  dplyr::select(c("sample", "total.reads"))

viral_len <- read.table("virus_len.txt", header = F)
colnames(viral_len) <- c("virus", "len")

# 3. join the total reads to each row based on sample
counts_long_plusreads <- left_join(counts_long,  total_reads_table, by = "sample")
# 4. join the gene length to each row based on gene
counts_long_plusreads_plustlen <- left_join(counts_long_plusreads, viral_len, by = "virus")

counts_long_plusreads_plustlen <- counts_long_plusreads_plustlen %>%
  dplyr::distinct()

# 5. compute the RPKM calculation and create count table from it 
rpkm_counts_final <- counts_long_plusreads_plustlen %>%
  dplyr::mutate(RPKM = (count*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("count", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) %>%
  column_to_rownames(var = "virus") %>%
  drop_na() %>%
  rownames_to_column(var = "Virus")
#######       #######      #######     ########     ######   ########    ########


# UPSET PLOT

# 1. Create presence abscence table
counts_pab <- rpkm_counts_final
counts_pab <- counts_pab %>%
  column_to_rownames(var = "Virus")
counts_pab[counts_pab > 0] <- 1

# () Get summary statistics on the number of vOTUs (for publication)
summary(colSums(counts_pab))
sd(colSums(counts_pab))

# 2. Get data ready for plot
# metadata ready
metadata_site<- metadata %>%
  dplyr::select(c("sample", "site")) 
#metadata_timepoints$sample <- gsub(metadata_timepoints$sample, pattern = "_", replacement = "-")

# presence abscence ready
pab_site <- counts_pab %>%
  rownames_to_column(var = "Virus") %>%
  pivot_longer(!Virus,names_to = "sample", values_to = "pab") %>%
  left_join(metadata_site, by = "sample") %>%
  dplyr::select(-c("sample")) %>%
  group_by(site, Virus) %>%
  dplyr::summarize(pab_total = sum(pab)) %>%
  pivot_wider(names_from = "site", values_from = "pab_total") %>%
  column_to_rownames(var = "Virus")
# make it to 0 and 1
pab_site[pab_site> 0] <- 1 

# UPSET PLOT TO SEE HOW MANY DNA EXPRESSION THEY SHARE
# https://github.com/hms-dbmi/UpSetR
library(UpSetR)
pdf(file="figures/virus_upset.pdf", width = 9, height = 5) # Change the path to your desired one to save it
upset(pab_site, sets = colnames(pab_site), order.by = "freq", empty.intersections = "on")
dev.off()

