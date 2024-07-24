# Analysis and visualization of KEGG annotations and functional traits

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



##############################################################################
######## SETUP
##############################################################################
# Set working directory 
setwd("Github/urban_greenspaces_aridcities/bacteria/data/")

# Loa data: 
ko_table <- read.table("ko_table.txt", header = T, row.names = 1)
ko_description <- read.table("ko_description.txt",  header = T, sep="\t")
ko_level1 <- read.table("ko_table_level1.txt", header = T, sep="\t", row.names = 1)
ko_level2 <- read.table("ko_table_level2.txt", header = T, sep="\t", row.names = 1)
ko_level3 <- read.table("ko_table_level3.txt", header = T, sep="\t", row.names = 1)
unknown <- read.table("ko_unknown.txt", header = T)
ko2level <- read.table("ko2level_Jan2021.txt", header = T, sep="\t" )
metadata <- read.csv("metadata.csv", row.names = 1)


#first let's get rid of the X at the beginning of the sample names
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

destroyX(ko_table)
destroyX(ko_level1)
destroyX(ko_level2)
destroyX(ko_level3)

# stablish rownames of unknown table as sample names
rownames(unknown) <- colnames(ko_table)


#sort all based on metadat rows
unknown <- unknown[rownames(metadata),]
ko_table <- ko_table[, rownames(metadata)]
ko_level1 <- ko_level1[, rownames(metadata)]
ko_level2 <- ko_level2[, rownames(metadata)]
ko_level3 <- ko_level3[, rownames(metadata)]

#check!!! 
all(rownames(metadata) == colnames(ko_table))
all(rownames(metadata) == colnames(ko_level1))
all(rownames(metadata) == colnames(ko_level2))
all(rownames(metadata) == colnames(ko_level3))
all(rownames(metadata) == rownames(unknown))

#record total KO counts per sample
metadata$total_abundance <- colSums(ko_table)




##############################################################################
######## GENES OF UNKNOWN FUNCTION
##############################################################################
#record proportion of unknown counts out of total counts per sample
metadata$prop.unknown<- unknown$unannotated_number/unknown$gene_number
metadata$prop.annotated.unknown<- unknown$annotated_unknown_number/unknown$gene_number
metadata$prop.annotated.general<- unknown$annotated_general_count/unknown$gene_number

anova(lm(prop.unknown~type, data= metadata))
anova(lm(prop.annotated.unknown~type, data= metadata))
summary(lm(prop.annotated.unknown~type, data= metadata))
anova(lm(prop.annotated.general~type, data= metadata))
summary(lm(prop.annotated.general~type, data= metadata))

ggplot(metadata, aes(y = prop.unknown,x = urban.natural, fill = urban.natural, color = urban.natural))+
  geom_boxplot(alpha=0.6)+
  geom_jitter()+ 
  xlab(NULL) +
  ylab("Proportion of unannanotated genes (%)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values = pal1)+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))


# STATS
Anova(lmer(metadata$prop.unknown ~ urban.natural + (1|site) , data = metadata), type=3)


##############################################################################
######## AVERAGE COPY NUMBER (acn)
##############################################################################

acn <- read.table("sample_16S_coverage.txt", header = T)
metadata$acn <- acn$V1

ggplot(data = metadata, aes(y = acn, x = urban.natural))+
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(fill= urban.natural, color = urban.natural))+
  geom_jitter(aes(color= urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average copy number")+
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_bw()+
  theme(legend.position = "none",  text = element_text(size=14))


# STATS
Anova(lmer(metadata$acn ~ urban.natural + (1|site) , data = metadata), type=3)



##############################################################################
######## AVERAGE GENOME SIZE (ags)
##############################################################################

ags <- read.table("ags_result.txt", header = T)
metadata$ags <- ags$ags

ggplot(data = metadata, aes(y = ags, x = urban.natural))+
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(fill= urban.natural, color = urban.natural))+
  #geom_boxplot(alpha=0.6, outlier.shape = NA)+
  geom_jitter(aes(color= urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average genome size")+
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_bw()+
  theme(legend.position = "none",  text = element_text(size=14))


# STATS
Anova(lmer(metadata$ags ~ metadata$urban.natural + (1|site), data = metadata), type = 3)



##############################################################################
######## GC CONTENT
##############################################################################
gc_mean <- read.table("gc_mean.txt")
gc_var <- read.table("gc_var.txt")

row.names(gc_mean)<- row.names(growth)
row.names(gc_var)<- row.names(growth)

#GC mean
ggplot(data = gc_mean, aes(y = V1, x = metadata$urban.natural))+
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(fill= metadata$urban.natural, color = metadata$urban.natural))+
  geom_jitter(aes(color = metadata$urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  ylab("Mean GC content (%)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))

metadata$gc_content <- gc_mean
anova(lmer(metadata$gc_content$V1 ~ urban.natural + (1|site) , data = metadata))
Anova(lmer(metadata$gc_content$V1 ~ urban.natural + (1|site) , data = metadata))
Anova(lmer(metadata$gc_content$V1 ~ urban.natural + (1|site) , data = metadata), type=3)
# likelihood ratio que compara el modelo con y sin random effect

#GC variance
ggplot(data = gc_var, aes(y = V1, x = metadata$urban.natural))+
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(fill= metadata$urban.natural, color = metadata$urban.natural))+
  geom_jitter(aes(color = metadata$urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  ylab("Variance GC content (%)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=14))

# STATS
metadata$gc_var <- gc_var
Anova(lmer(metadata$gc_var$V1 ~ metadata$urban.natural + (1|site), data = metadata), type = 3)




##############################################################################
######## CODON USAGE BIAS
##############################################################################
growth <- read.table("growth_coverage.txt", header = T)

# Mean of the codon usage bias of each highly expressed gene relative to all other genes
metadata$cubhe <- growth$cubhe
ggplot(data = growth, aes(y = cubhe, x = metadata$urban.natural))+
  geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(fill= metadata$urban.natural, color = metadata$urban.natural))+
  geom_jitter(aes(color =  metadata$urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average codon usage bias")+
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("codonusagebias.pdf", device = "pdf", width = 4, height = 4 , units = "in")

Anova(lmer(metadata$cubhe ~ urban.natural + (1|site) , data = metadata), type=3)








##############################################################################
######## SUGAR-ACID metabolism
##############################################################################

# calculate abundance of genes in sugar vs. acid metabolism and add to metadata
sugar_acid <- read.csv("sugar_acid_ko_gralka23.csv")

# list of ribosomal kos
scgs <- c("K02950", "K02992", "K02906", "K02926", "K02890", "K02874", "K02931", "K02994", "K02952", "K02948", "K02871", "K02961", "K02881", "K02986")

# obtain the counts of the ribosomal kos
scgs_counts <- ko_table %>%
  rownames_to_column(var = "ko") %>%
  filter(ko %in% scgs) %>%
  column_to_rownames(var = "ko")

# calculate average counts
av_scgs <- as.data.frame(colMeans(scgs_counts))

av_scgs <- av_scgs %>%
  rownames_to_column(var = "sample") %>%
  mutate(sample = gsub(sample, pattern = "^X", replacement = ""))
colnames(av_scgs) <- c("sample", "av_counts_ribosomal")


# normalize all ko counts by the av_counts of these ribosomal genes
ko_table_norm <- ko_table %>%
  rownames_to_column(var = "ko") %>% 
  pivot_longer(-ko, names_to = "sample", values_to = "counts") %>%
  mutate(sample = gsub(sample, pattern = "^X", replacement = "")) %>%
  left_join(., av_scgs, by = "sample") %>%
  mutate(norm_counts = counts/av_counts_ribosomal) %>%
  dplyr::select(-c("counts", "av_counts_ribosomal")) %>%
  pivot_wider(names_from = "sample", values_from = "norm_counts") %>%
  column_to_rownames(var = "ko")


sugar_acid_ko_count <- ko_table_norm %>%
  rownames_to_column(var = "ko") %>%
  left_join(., sugar_acid, by = "ko") %>%
  na.omit() %>%
  dplyr::select(-c("metabolism", "path", "description")) %>%
  dplyr::select(-ko) %>%
  group_by(source) %>%
  summarize_all(sum) 

sugar_acid_ko_count_to_join <- as.data.frame(t(sugar_acid_ko_count))
colnames(sugar_acid_ko_count_to_join) <- sugar_acid_ko_count_to_join[1,] 
sugar_acid_ko_count_to_join <- sugar_acid_ko_count_to_join[-1,]

sugar_acid_ko_count_to_join <- sugar_acid_ko_count_to_join %>%
  rownames_to_column(var = "sample")

metadata <- left_join(metadata,sugar_acid_ko_count_to_join,by = "sample" )


# calculate the ratio sugar/acid
metadata$sugar_acid_ratio <- as.numeric(metadata$sugar)/as.numeric(metadata$acid)

# boxplot
ggplot(data = metadata, aes(y = sugar_acid_ratio , x = urban.natural))+
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(fill= urban.natural, color = urban.natural))+
  geom_jitter(aes(color= urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Ratio of abundance genes in the metabolism of sugars vs acids")+
  scale_x_discrete(labels=c('Natural', 'Urban greenspaces'))+
  theme_bw()+
  theme(legend.position = "none",  text = element_text(size=16))


# STATS
Anova(lmer(metadata$sugar_acid_ratio ~ urban.natural + (1|site) , data = metadata), type=3)


library(SparkR)
a = -2.023206
s =  6.080268
metadata$SAP <- tanh(s*(as.numeric(metadata$sugar)/metadata$total.reads) + a*(as.numeric(metadata$acid)/metadata$total.reads))

metadata_urban <- metadata %>%
  dplyr::filter(landuse == "urban")

metadata_natural<- metadata %>%
  dplyr::filter(!landuse == "urban")

sum(as.numeric(metadata_urban$sugar))

urban_sap <- tanh(s*sum(as.numeric(metadata_urban$sugar)) + a*sum(as.numeric(metadata_urban$acid)))
natural_sap <- tanh(s*sum(as.numeric(metadata_natural$sugar)) + a*sum(as.numeric(metadata_natural$acid)))

# boxplot
ggplot(data = metadata, aes(y = SAP, x = urban.natural))+
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(fill= urban.natural, color = urban.natural))+
  geom_jitter(aes(color= urban.natural), size = 2.5)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Sugar-acid preference (SAP)")+
  scale_x_discrete(labels=c('Natural', 'Urban'))+
  theme_bw()+
  theme(legend.position = "none",  text = element_text(size=16))

# STATS
Anova(lmer(metadata$SAP ~ urban.natural + (1|site) , data = metadata), type=3)




