# Soil physicochemical properties analysis


# Load packages
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(metagenomeSeq) #bioconductor
library(gRodon)
library(Biostrings)
library(ggrepel)
library(lme4)
library(tidyverse)
library(ggbiplot)
library(car)
library(performance)
pal1 = c("#FFC20A", "#0C7BDC")




## Soil chemistry

# Set working directory 
setwd("Github/urban_greenspaces_aridcities/bacteria/data/")

chem <- read.table('soil_chemistry.txt', header = T)
metadata <- read.csv("metadata.csv", row.names = 1)

chem <- chem %>%
  arrange(Lab)

chem$Lab <- gsub(chem$Lab, pattern = "-T2-", replacement = ".T1.")
chem$Lab <- gsub(chem$Lab, pattern = "M", replacement = "")

chem <- chem %>%
  dplyr::select(-c("SampleID", "Date", "Crop", "FreeLime", "id")) %>%
  column_to_rownames(var = "Lab")

#do pca and visualize
chempca <- prcomp(chem, center = TRUE, scale = TRUE)

ggbiplot(chempca, groups = metadata$urban.natural, ellipse = T)+
  #geom_text(label = chem$id, size = 2)+
  geom_point(aes(colour=metadata$urban.natural), size = 4) +
  scale_color_manual(values=pal1)+
  #scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_bw() +
  theme(legend.title = element_blank())

metadata <- metadata %>%
  mutate(color = case_when(urban.natural == "urban" ~ "#0C7BDC", 
                           TRUE ~ "#FFC20A"))


# Z normalize
chem_scaled <- as.data.frame(scale(chem))
## Pheatmap
pheatmap(as.matrix(chem_scaled), 
         display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 8, 
         cluster_cols = FALSE, 
         color = hcl.colors(50, "BluYl"))



