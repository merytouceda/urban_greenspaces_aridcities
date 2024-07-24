# phatype output parsing for land use project

library(tidyverse)
library(dplyr)
library(lme4) #mixed effect models
library(lmerTest)
library(agricolae)
library(car)
library(performance)
pal1 = c("#FFC20A", "#0C7BDC")
pal2 = c("#009E73","#E66100", "#5D3A9B")

setwd("Github/urban_greenspaces_aridcities/virus/data")
phatype <- read.csv("landuse_phatype_out.csv", header = T)
metadata<- read.csv("metadata.csv", row.names = 1)
counts <- read.table("viral_species_count_table.txt", header = T, row.names = 1, sep = "\t")
# remove X from before sample name
colnames(counts) <- gsub(colnames(counts), pattern = "^X", replacement = "")
colnames(counts) <- gsub(colnames(counts), pattern = ".Read.Count", replacement = "")


head(phatype)
summary(phatype$Score)

# count cases of virulent and temperate based on landuse 
phatype_summary <- phatype %>%
  filter(Score > 0.8) %>%
  mutate(
    Landuse = case_when(
      str_detect(Contig, "^HP") | str_detect(Contig, "^RP") ~ "Urban", 
      TRUE ~ "Natural")) %>%
  group_by(Landuse, Pred) %>%
  tally() # like count but you have to do the grouping


# count by sample
phatype_sample<- phatype %>%
  separate(Contig, c("Sample", "Contig"), sep= "_") %>%
  filter(Score > 0.8) %>%
  group_by(Sample, Pred) %>%
  tally() %>%
  mutate(
    Landuse = case_when(
      str_detect(Sample, "^HP-") | str_detect(Sample, "^RP") ~ "Urban", 
      TRUE ~ "Natural")) %>%
  mutate(
    Vegetation = case_when(
      str_detect(Sample, "^SC") ~ "Shrubland", 
      str_detect(Sample, "^RC") ~ "Forest",
      TRUE ~ "Grassland"))


# LIFE STYLE PREVALENCE
# boxplot with facet wrap (could do this with pointrange instead) 
ggplot(phatype_sample, aes(x = Landuse, y = n))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(color = Landuse, fill = Landuse))+ 
  geom_jitter(aes(color= Landuse), size = 2.5)+
  xlab(NULL) +
  ylab("Prevalence (number of viral genomes)")+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))+
  facet_wrap(~Pred)


# boxplot for supplementary
ggplot(phatype_sample, aes(x = Landuse, y = n))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA)+ 
  geom_jitter(aes(color= Vegetation), size = 2.5)+
  xlab(NULL) +
  ylab("Prevalence (number of viral genomes)")+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))+
  facet_wrap(~Pred)


# stats
# t.test based on landuse
t.test(phatype_sample$n ~ phatype_sample$Landuse,subset = phatype_sample$Pred == "virulent")
t.test(phatype_sample$n ~ phatype_sample$Landuse,subset = phatype_sample$Pred == "temperate")

# t.test based on landuse
t.test(phatype_sample$n ~ phatype_sample$Pred, subset = phatype_sample$Landuse == "Urban")
t.test(phatype_sample$n ~ phatype_sample$Pred, subset = phatype_sample$Landuse == "Natural")


# anova
anova(lm(phatype_sample$n ~ phatype_sample$Landuse,subset = phatype_sample$Pred == "virulent"))
anova(lm(phatype_sample$n ~ phatype_sample$Landuse,subset = phatype_sample$Pred == "temperate"))

anova(lm(phatype_sample$n ~ phatype_sample$Pred, subset = phatype_sample$Landuse == "Urban"))
anova(lm(phatype_sample$n ~ phatype_sample$Pred, subset = phatype_sample$Landuse == "Natural"))

# anova
anova(lm(phatype_sample$n ~ phatype_sample$Vegetation,subset = phatype_sample$Pred == "virulent"))
anova(lm(phatype_sample$n ~ phatype_sample$Vegetation,subset = phatype_sample$Pred == "temperate"))

anova(lm(phatype_sample$n ~ phatype_sample$Vegetation, subset = phatype_sample$Landuse == "Urban"))
anova(lm(phatype_sample$n ~ phatype_sample$Vegetation, subset = phatype_sample$Landuse == "Natural"))




############################################################
# make counts into RPKM normalized
############################################################
# load length
viral_len <- read.table("virus_len.txt", header = F)
colnames(viral_len) <- c("Contig", "Length")
viral_len <- unique(viral_len)

# create total reads
#total_reads <- as.data.frame(colnames(counts))
#total_reads$total.reads <- colSums(counts)
#colnames(total_reads) <- c("Sample", "total.reads")

total_reads_table <- metadata %>%
  dplyr::select(c("sample", "total_abundance"))
colnames(total_reads_table) <- c("Sample", "total.reads")

counts_long <- counts %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(cols = -c("Contig"), values_to = "Count", names_to= "Sample")

# 3. join the total reads to each row based on sample
counts_long_plusreads <- left_join(counts_long, total_reads_table, by = "Sample")
# 4. join the gene length to each row based on gene
counts_long_plusreads_plustlen <- left_join(counts_long_plusreads, viral_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_counts <- counts_long_plusreads_plustlen  %>%
  dplyr:: mutate(RPKM = (Count*10e6)/(total.reads*as.numeric(Length))) %>%
  dplyr::select(-c("Count", "total.reads", "Length")) %>%
  pivot_wider(names_from = Sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig")



phatype_counts <- rpkm_counts %>%
  rownames_to_column(var = "Contig") %>%
  left_join(phatype, by = "Contig") %>%
  dplyr::select(-Score) %>%
  na.omit() %>% 
  dplyr::select(-Contig) %>%
  pivot_longer(-Pred, names_to = "Sample", values_to = "Count") %>%
  dplyr::group_by(Sample, Pred) %>%
  dplyr::summarize(n = sum(Count)) %>%
  dplyr::mutate(
    landuse = case_when(
      str_detect(Sample, "^HP") | str_detect(Sample, "^RP") ~ "Urban", 
      TRUE ~ "Natural")) %>%
  dplyr::mutate(
    Vegetation = case_when(
      str_detect(Sample, "^SC") ~ "Shrubland", 
      str_detect(Sample, "^RC") ~ "Forest",
      TRUE ~ "Grassland"))

# by landuse
ggplot(phatype_counts, aes(x = landuse, y = n))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(color = landuse, fill = landuse))+ 
  geom_jitter(aes(color= landuse), size = 2)+
  xlab(NULL) +
  ylab("Abundance in RPKM")+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))+
  facet_wrap(~Pred)


phatype_counts <- phatype_counts %>%
  dplyr::mutate(site = str_extract(Sample, pattern = "[^.]*"))


phatype_counts_u <- phatype_counts %>%
  dplyr::filter(landuse == "Urban")
phatype_counts_n <- phatype_counts %>%
  dplyr::filter(!landuse == "Urban")

anova(lmer(phatype_counts_n$n ~ Pred + (1|site) , data = phatype_counts_n))
anova(lmer(phatype_counts_n$n ~ Pred + (1|site) , data = phatype_counts_n), type = 3)

anova(lmer(phatype_counts_u$n ~ Pred + (1|site) , data = phatype_counts_u))
anova(lmer(phatype_counts_u$n ~ Pred + (1|site) , data = phatype_counts_u), type = 3)

phatype_counts_t <- phatype_counts %>%
  dplyr::filter(Pred == "temperate")

phatype_counts_v <- phatype_counts %>%
  dplyr::filter(!Pred == "temperate")

anova(lmer(phatype_counts_t$n ~ landuse + (1|site) , data = phatype_counts_t))
anova(lmer(phatype_counts_t$n ~ landuse + (1|site) , data = phatype_counts_t), type = 3)
anova(lmer(phatype_counts_v$n ~ landuse + (1|site) , data = phatype_counts_v), type = 3)

# by habitat
ggplot(phatype_counts, aes(x = landuse, y = n))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA)+ 
  geom_jitter(aes(color= Vegetation), size = 2)+
  xlab(NULL) +
  ylab("Abundance in RPKM")+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))+
  facet_wrap(~Pred)




# Proportion of life strategy annotated viruses: 1447/4306 = 0.336
############################################################
# lytic/lysogenic ratio 
############################################################
# natural: 
nat_summary_phatype_counts <- phatype_counts %>%
  dplyr::filter(!landuse == "Urban") %>%
  dplyr::group_by(Pred) %>%
  dplyr::summarize(average_abundance = mean(n))
# 7.820301/14.416486 = 0.54

# urban: 
urb_summary_phatype_counts <- phatype_counts %>%
  dplyr::filter(landuse == "Urban") %>%
  dplyr::group_by(Pred) %>%
  dplyr::summarize(average_abundance = mean(n))
# 20.70406/30.56493 = 0.6773796

############################################################
# Abundance in general of viruses
############################################################
colSums(rpkm_counts)

metadata$total_viral_abundance <- colSums(rpkm_counts)
ggplot(metadata, aes(x = urban.natural, y = total_viral_abundance))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA, aes(color = urban.natural, fill = urban.natural))+ 
  geom_jitter(aes(color= urban.natural), size = 2)+
  xlab(NULL) +
  ylab("Abundance in RPKM")+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))


anova(lmer(total_viral_abundance ~ urban.natural + (1|site) , data = metadata))


ggplot(metadata, aes(x = urban.natural, y = total_viral_abundance))+ 
  geom_boxplot(alpha=0.6, outlier.shape = NA)+ 
  geom_jitter(aes(color= vegetation_structure), size = 2)+
  xlab(NULL) +
  ylab("Abundance in RPKM")+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))


anova(lm(total_viral_abundance ~ vegetation_structure + (1|site) , data = metadata))

############################################################
# % of viral reads out of all reads
############################################################


total_reads <- metadata %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "total.reads"))

sum(colSums(counts))/sum(total_reads$total.reads)

# 0.006429418 = 6.4%

