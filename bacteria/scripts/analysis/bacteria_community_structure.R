# Bacterial community structure visualization and analysis

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

setwd("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/bacteria/data")
bracken_lu <- read.csv("braken_all_samples.csv", header = T, row.names = 2)
bracken_lu$X <- NULL
bracken_lu <- as.data.frame(t(bracken_lu))
metadata <- read.csv("../../metadata.csv", row.names = 1)

# Make sure  that the row names are in the same order as in metadata
rownames(bracken_lu) <- gsub(rownames(bracken_lu), pattern = "^X", replacement = "")
rownames(bracken_lu) == rownames(metadata)

# Replace nas for 0s
bracken_lu <- bracken_lu %>% replace(is.na(.), 0)

# Filter low prevalence and abundance cases

# 1. Filter based on PREVALENCE (keep those present in > 10% of samples)
# create a presence abscence table
bracken_lu_pab <- as.data.frame(t(bracken_lu))
bracken_lu_pab[bracken_lu_pab > 0] <- 1

# Keep those columns (species) that are prevalent
# meaning, their abundance has to be > 0 in at least 10% of the samples
a <- dim(bracken_lu)[1]*0.1 # 10% of samples
prev_bracken <- bracken_lu_pab[,colSums(bracken_lu_pab) >= a] # take out those that don't have a 1 in more than 10%
bracken_lu <- bracken_lu[,colnames(bracken_lu) %in% colnames(prev_bracken)] # match the original table to the one with the prevalent species


# 2. Filter based on ABUNDANCE (keep those bacteria whose total abundance >= 0.00001)
ab_bracken_lu <- bracken_lu[,colSums(bracken_lu)/sum(colSums(bracken_lu)) >= 0.00001]



####################################################################################
############## GENERAL ECOLOGY
####################################################################################

# ALPHA DIVERSITY
summary(rowSums(ab_bracken_lu))
metadata$bact.species.rich.rar <- specnumber(rrarefy(ab_bracken_lu, sample = min(rowSums(ab_bracken_lu))))

# boxplot
ggplot(metadata, aes(x = urban.natural, y = bact.species.rich.rar))+
  geom_boxplot(alpha=0.6, outlier.shape = NA)+
  geom_jitter(aes(color= vegetation_structure), size = 3)+
  xlab(NULL) + 
  ylab("Number of different viral species")+
  scale_fill_manual(values=pal1)+
  #scale_color_manual(values=pal2)+
  #scale_fill_see()+
  scale_color_manual(values = pal2)+
  #scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/bacteria/figures/bact_rich.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# STATS
Anova(lmer(bact.species.rich.rar ~ urban.natural*landuse2 + (1|site) , data = metadata), type = 3)


# BETA-DIVERSITY
# Calculate bray curtis dissimilarity
counts.bray <- vegdist(ab_bracken_lu, method="bray")

#Generate ordination with NMDS (non-multidimensional scaling)
counts.nmds <- metaMDS(counts.bray, k=2, try = 100)
metadata$Axis01 = counts.nmds$points[,1]
metadata$Axis02 = counts.nmds$points[,2]
counts.nmds$stress #0.094


# The only option is to change the shape ATR and the color for plot! 
ggplot(metadata, aes(Axis01, Axis02))+
  geom_point(aes(color= vegetation_structure), size=3)+
  stat_ellipse(aes(group = urban.natural), linetype = 2)+
  #geom_text(aes(label = metadata$sample))+
  #geom_text_repel(aes(label = Group, color = Category, alpha = 0.5))+
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Category, color = Category))+
  scale_color_manual(values=pal2)+
  theme_classic()+
  theme(legend.position="none", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/GitHub/urban_greenspaces_aridcities/bacteria/figures/bact_commcomp.pdf", device = "pdf", width = 5, height = 5 , units = "in")

# STATS
adonis2(counts.bray ~ urban.natural*landuse2, data = metadata, permutations = 999, method = "bray")


# beta dispersion
mod <- betadisper(counts.bray, as.factor(metadata$urban.natural))
# mod plot
plot(mod)
# boxplot
boxplot(mod)
# STATS





# calculate average per sample species
#bacteria
ab_bracken_lu_pab <- as.data.frame(t(ab_bracken_lu))
ab_bracken_lu_pab[ab_bracken_lu_pab > 0] <- 1

mean(rowSums(ab_bracken_lu_pab))
sd(rowSums(ab_bracken_lu_pab))


#UPSET PLOT

metadata_site<- metadata %>%
  dplyr::select(c("sample", "site")) 
#metadata_timepoints$sample <- gsub(metadata_timepoints$sample, pattern = "_", replacement = "-")

pab_site <- ab_bracken_lu_pab %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(!sample,names_to = "bacteria", values_to = "pab") %>%
  left_join(metadata_site, by = "sample") %>%
  dplyr::select(-c("sample")) %>%
  group_by(site, bacteria) %>%
  dplyr::summarize(pab_total = sum(pab)) %>%
  pivot_wider(names_from = "site", values_from = "pab_total") %>%
  column_to_rownames(var = "bacteria")

# make it to 0 and 1
pab_site[pab_site> 0] <- 1 

# UPSET PLOT TO SEE HOW MANY DNA EXPRESSION THEY SHARE
# https://github.com/hms-dbmi/UpSetR
library(UpSetR)
pdf(file="/figures/bact_upset.pdf", width = 9, height = 5) # Change the path to your desired one to save it
upset(pab_site, sets = colnames(pab_site), order.by = "freq", empty.intersections = "on")
dev.off()