# kegg analysis of plasmids and mges from urban project


library(tidyverse)
library(see)
library(pals)
library(vegan)
library(Polychrome)
library(car)
library(lme4)
library(ggExtra)
library(lmPerm)
library (performance)
pal1 <- glasbey(28)

pal2 = c("#56282D", "#77966D")

#ko2level <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/ko2level_Jan2021.csv")
ko2level <- read_csv("/Volumes/BunnyBike/mge_urban/local/kegg_ko_complete.csv")

# load table of kegg annotation on plasmids
plasmids_kegg <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plasmids_kegg_results.csv")

# add the level information to the annotated table
plasmids_kegg_level <- plasmids_kegg %>%
  mutate(KO = query_name) %>%
  left_join(., ko2level, by = "KO")

# keep only the best hit to each plasmid gene
plasmids_kegg_level_best <- plasmids_kegg_level %>%
  dplyr::group_by(target_name) %>% # group by same read name
  dplyr::arrange(e_value,.by_group = T) %>% # sort by evalue within group
  slice_head(n = 1) # select the first hit (best one) for each group 


# see how many lack category
not_categorized <-plasmids_kegg_level %>% # should be as little as possible, otherwise the ko2level might not be updated
  filter(is.na(level1))


# how many xenobiotic and what xenobiotics
xargs <- plasmids_kegg_level_best %>%
  filter(level2 %in% c("09111 Xenobiotics biodegradation and metabolism")) %>%
  dplyr::group_by(level3) %>%
  dplyr::summarize(n = n())

ggplot(xargs, aes(x = reorder(level3, n), y = n))+
  geom_col()+
  xlab("")+
  ylab("Number of genes")+
  coord_flip()
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/xargs_compounds.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# by soil environment
ptu_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_count_table.txt", sep = "\t", header = T)

ptu_pab <- ptu_counts %>%
  column_to_rownames(var = "Contig") %>%
  rownames_to_column(var = "plasmid") %>%
  pivot_longer(-c("plasmid"), names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(pab = case_when(count ==  0 ~ 0, 
                         TRUE ~ 1)) %>%
  select(-c("count")) %>%
  pivot_wider(names_from = "sample", values_from = "pab")

xargs_tojoin <- plasmids_kegg_level_best %>%
  filter(level2 %in% c("09111 Xenobiotics biodegradation and metabolism")) %>%
  ungroup() %>%
  mutate(plasmid = gsub(target_name, pattern = "_[0-9]*$", replacement = "")) %>%
  select(c("plasmid", "level3"))

xargs_counts_class <-xargs_tojoin  %>%
  dplyr::left_join(., ptu_pab, by = "plasmid") %>%
  drop_na() %>%
  dplyr::select(-plasmid) %>%
  pivot_longer(-level3, names_to = "sample", values_to="count") %>%
  dplyr::group_by(level3, sample) %>%
  dplyr::summarize(counts = sum(count)) %>%
  dplyr::mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural"))

ggplot(xargs_counts_class, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~level3)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number of genes")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/xargs_class_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")





##########
# See abundance of different functions
##########
# convert variables of interest into factors
levels(as.factor(plasmids_kegg_level$level1))
levels(as.factor(plasmids_kegg_level$level2))
levels(as.factor(plasmids_kegg_level$level3))

# Make a table with the overall counts for all level 2 categories
count_level2 <- plasmids_kegg_level_best %>%
  dplyr::select(c("KO", "target_name", "level2")) %>%
  dplyr::mutate(plasmid = gsub(target_name, pattern = "_[0-9]*$", replacement = "")) %>%
  dplyr::mutate(count = 1) %>%
  dplyr::group_by(plasmid, level2) %>%
  dplyr::summarize(counts = sum(count)) %>%
  dplyr::filter(!level2 %in% c("09167 Endocrine and metabolic disease", "09166 Cardiovascular disease", "09165 Substance dependence",
                        "09161 Cancer: overview", "09162 Cancer: specific types", "09163 Immune disease", 
                        "09164 Neurodegenerative disease", "09149 Aging","09151 Immune system" , "09152 Endocrine system" ,
                        "09153 Circulatory system","09154 Digestive system" ,"09155 Excretory system" ,"09156 Nervous system",
                        "09157 Sensory system", "09158 Development and regeneration", "09144 Cellular community - eukaryotes" )) %>%
  dplyr::mutate(level2 = gsub(level2, pattern = "[0-9]*", replacement = ""))  %>%
  dplyr::mutate(level2 = case_when(str_detect(level2, "Unclassified") ~ "Unclassified", 
                            TRUE ~ level2))

# make percenteges table to see relative abundance
count_level2_percentages <- count_level2 %>%
  group_by(plasmid) %>%
  mutate(percentage = (counts / sum(counts)) * 100) %>%
  ungroup() %>%
  mutate(percentage = round(percentage, 2))


# Barplot counts
ggplot(count_level2_percentages, aes ( x = plasmid, y = counts, fill = level2)) + 
  geom_col()+
  scale_fill_manual(values = pal1)+
  xlab("Plasmid")+
  ylab("Number of genes")+
  theme(legend.position = "bottom", axis.ticks.x=element_blank(), axis.text.x=element_blank())
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmids_functions_kegg.pdf", device = "pdf", width = 15, height = 7 , units = "in")


# Barplot %s
ggplot(count_level2_percentages, aes ( x = plasmid, y = percentage, fill = level2)) + 
  geom_col()+
  scale_fill_manual(values = pal1)+
  xlab("Plasmid")+
  ylab("Number of genes")+
  theme(legend.position = "bottom", axis.ticks.x=element_blank(), axis.text.x=element_blank())



##################
# Add information about plasmids (conjugation genes, taxonomy)
##################

all_summary_filtered <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/circular_10kb_plasmids_summary.csv")

all_summary_filtered_join <- all_summary_filtered %>%
  select(c("seq_name", "conjugation_genes", "amr_genes"))

# add it to kegg counts
plasmids_kegg_level_best_info <- plasmids_kegg_level_best %>%
  mutate(seq_name = gsub(target_name, pattern = "_[0-9]*$", replacement = "")) %>%
  left_join(., all_summary_filtered_join, by = "seq_name")

# Prep for barplot but only conjugative plasmids
count_level2_conj <- plasmids_kegg_level_best_info %>%
  dplyr::filter(!is.na(conjugation_genes)) %>%
  dplyr::select(c("KO", "target_name", "level2")) %>%
  dplyr::mutate(plasmid = gsub(target_name, pattern = "_[0-9]*$", replacement = "")) %>%
  dplyr::mutate(count = 1) %>%
  dplyr::group_by(plasmid, level2) %>%
  dplyr::summarize(counts = sum(count)) %>%
  dplyr::filter(!level2 %in% c("09167 Endocrine and metabolic disease", "09166 Cardiovascular disease", "09165 Substance dependence",
                        "09161 Cancer: overview", "09162 Cancer: specific types", "09163 Immune disease", 
                        "09164 Neurodegenerative disease", "09149 Aging","09151 Immune system" , "09152 Endocrine system" ,
                        "09153 Circulatory system","09154 Digestive system" ,"09155 Excretory system" ,"09156 Nervous system",
                        "09157 Sensory system", "09158 Development and regeneration", "09144 Cellular community - eukaryotes" )) %>%
  dplyr::mutate(level2 = gsub(level2, pattern = "[0-9]*", replacement = ""))  %>%
  dplyr::mutate(level2 = case_when(str_detect(level2, "Unclassified") ~ "Unclassified", 
                            TRUE ~ level2))
# Barplot
ggplot(count_level2_conj, aes ( x = plasmid, y = counts, fill = level2)) + 
  geom_col()+
  scale_fill_manual(values = pal1)+
  xlab("Plasmid")+
  ylab("Number of genes")+
  theme(legend.position = "bottom", axis.ticks.x=element_blank(), axis.text.x=element_blank())
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmids_functions_kegg_conj.pdf", device = "pdf", width = 15, height = 5 , units = "in")


# Barplot of level2 functions not separated by plasmids
level_2_number <- count_level2_conj %>%
  group_by(level2) %>%
  summarize(Counts = sum(counts)) %>%
  drop

ggplot(level_2_number, aes(x = reorder(level2, Counts), y = Counts)) + 
  geom_col()+
  xlab("")+
  ylab("Number of annotated genes")+
  coord_flip()+
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/functions_kegg.pdf", device = "pdf", width = 7, height = 5 , units = "in")






##################
# Counts
##################

############################################ PTU counts
gene_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plasmid_genes_count_table.txt", sep = "\t", header = T)
ptu_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_count_table.txt", sep = "\t", header = T)

md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv") %>%
  arrange(desc(sample))

kegg_counts<- plasmids_kegg_level_best_info %>%
  mutate(Contig = seq_name) %>%
  left_join(., ptu_counts, by = "Contig") %>%
  filter(!is.na(X11B.T1.3.Read.Count))
  
kegg_gene_counts <-  plasmids_kegg_level_best_info %>%
  mutate(Contig = target_name) %>%
  left_join(., gene_counts, by = "Contig") %>%
  filter(!is.na(X11B.T1.3.Read.Count))


# richness of genes
gene_counts_rich <- gene_counts %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  #mutate(sample = gsub(sample, pattern = "\\.", replacement = "-")) %>%
  pivot_wider(names_from = "Contig", values_from = "count") %>% 
  #arrange(desc(sample)) %>%
  column_to_rownames(var = "sample")

# make sure samples are in the same order
rownames(gene_counts_rich) == md$sample

# Caluclate richness
summary(rowSums(gene_counts_rich))
md$p_gene_rich <- specnumber(rrarefy(gene_counts_rich, sample = 4863))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "p_gene_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(p_gene_rich),
                   p_gene_rich = mean(p_gene_rich))

#pal1 = c("#FF823B", "#359951")

ggplot(md, aes(x = urban.natural, y = p_gene_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  p_gene_rich-sd, ymax =  p_gene_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("Plasmid functional richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_gene_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(p_gene_rich ~ urban.natural + (1|site), data = md), type = 3)

## See diagnostics at the end



# functional diversity / ptu diversity
md <- md %>%
  arrange(desc(sample))

md_final_join <- md_final %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "ptu_rich"))

md$sample == rownames(md_final)

md<- md %>%
  left_join(., md_final_join, by = "sample")
  mutate(ratio = p_gene_rich/ptu_rich) %>%
  
# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "ratio")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(ratio),
                   ratio = mean(ratio))

ggplot(md, aes(x = urban.natural, y = ratio)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  ratio-sd, ymax =  ratio+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Ratio between functional richness / plasmid richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_gene_ptu_richness_ratio.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(ratio ~ urban.natural + (1|site), data = md), type = 3)
#urban.natural   6.1142  1    0.01341 *

library(ggpubr)
ggplot(md, aes(x = ptu_rich, y = p_gene_rich))+
  geom_point(aes(color = urban.natural))+
  geom_smooth(method = lm, aes(color = urban.natural, fill = urban.natural))+
  stat_cor(method = "pearson", aes(color = urban.natural))+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab("PTU richness")+
  ylab("Functional richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 16))




#### Heatmap presence abscence
# Make a table that tells me if a function is present in urban, natural, or both 
ptu_pab <- ptu_counts %>%
  pivot_longer(-c("Contig"), names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",  "HP.T1.6",   "HP.T1.9",
                                                 "RC.T1.3",   "RC.T1.6",   "RC.T1.9" ) ~ "Urban", 
                                   TRUE ~ "Natural")) %>%
  group_by(Contig, urban.natural) %>%
  summarize(counts = sum(count)) %>%
  mutate(pab = case_when(counts == 0 ~ 0, 
                         TRUE ~ 1)) %>%
  group_by(Contig) %>%
  mutate(presence = case_when(
    sum(pab) == 2 ~ "Both",
    sum(pab == 1 & urban.natural == "Urban") == 1 ~ "Urban only",
    sum(pab == 1 & urban.natural == "Natural") == 1 ~ "Natural only"
  )) %>%
  ungroup() %>%
  distinct(Contig, presence) %>%
  drop_na()

# Add the presence abscence information to a level2 and contig table
kegg_pab <- kegg_counts %>%
  ungroup() %>%
  select(c("level2", "Contig")) %>%
  #mutate(Contig = gsub(Contig, pattern = "_[0-9]*$", replacement = "")) %>% # only use this if contig has the gene _1 at the end
  left_join(., ptu_pab, by = "Contig")%>%
  filter(!level2 %in% c("09167 Endocrine and metabolic disease", "09166 Cardiovascular disease", "09165 Substance dependence",
                        "09161 Cancer: overview", "09162 Cancer: specific types", "09163 Immune disease", 
                        "09164 Neurodegenerative disease", "09149 Aging","09151 Immune system" , "09152 Endocrine system" ,
                        "09153 Circulatory system","09154 Digestive system" ,"09155 Excretory system" ,"09156 Nervous system",
                        "09157 Sensory system", "09158 Development and regeneration", 
                        "09144 Cellular community - eukaryotes", "09194 Poorly characterized" )) %>%
  drop_na() #%>%
# this removes the order which is bad
 # mutate(level2 = gsub(level2, pattern = "^[0-9]* ", replacement = ""))


# Tile plot
ggplot(kegg_pab, aes(x = Contig, y = level2,fill = presence)) +
  geom_tile() +
  scale_fill_manual(values = c("Both" = "skyblue", 
                               "Urban only" = "#359951",
                               "Natural only" = "#FF823B")) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Plasmid", y = "", fill = "")
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_functions.pdf", device = "pdf", width = 11, height = 9 , units = "in")





### Heatmap with hirearchichal clustering
# Create a binary matrix for clustering
cluster_matrix <- kegg_pab %>%
  distinct(Contig, level2, presence) %>%
  mutate(presence_binary = case_when(
    presence == "Both" ~ 1,
    presence == "Urban only" ~ 2,
    presence == "Natural only" ~ 3
  )) %>%
  select(Contig, level2, presence_binary) %>%
  pivot_wider(names_from = level2, values_from = presence_binary, values_fill = 0) %>%
  column_to_rownames("Contig")

# Perform hierarchical clustering on contigs (columns in the plot)
contig_order <- hclust(dist(cluster_matrix))$order
contig_levels <- rownames(cluster_matrix)[contig_order]

# Optionally, also cluster the level2 categories (rows)
level2_order <- hclust(dist(t(cluster_matrix)))$order
level2_levels <- colnames(cluster_matrix)[level2_order]

# Apply the ordering to your plot
kegg_pab_ordered <- kegg_pab %>%
  mutate(Contig = factor(Contig, levels = contig_levels),
         level2 = factor(level2, levels = level2_levels))

# Plot with ordered axes
ggplot(kegg_pab_ordered, aes(x = Contig, y = level2, fill = presence)) +
  geom_tile() +
  scale_fill_manual(values = c("Both" = "skyblue", 
                               "Urban only" = "#359951",
                               "Natural only" = "#FF823B")) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Plasmid", y = "", fill = "")
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_functions_hclust.pdf", device = "pdf", width = 11, height = 9 , units = "in")



##########
## Richness
##########
levels(kegg_pab$Contig)

# see what Contig-KO pairs are present multiple times
df <- kegg_counts  %>%
  group_by(KO, Contig) %>%
  tally() %>%
  arrange(desc(n))

# remove duplicated pairs
kegg_counts_rich <- kegg_counts %>%
  ungroup() %>%
  select(-c("kegg_id" , "target_name" ,"target_accession" , "query_name" ,"query_accession","e_value","score" , "bias" , "best_domain_e_value" ,
            "best_domain_score" ,"best_domain_bias", "...1" , "level1" ,"level2" ,"level3" ,"KO_description" ,  "seq_name", "conjugation_genes",
            "amr_genes")) %>%
  dplyr::distinct(KO, Contig, .keep_all = TRUE) %>%
  select(-Contig) 

names_columnas <- c("KO","11B.T1.3" , "11B.T1.6" , "11B.T1.9" , 
                    "8.T1.3",    "8.T1.6" ,   "8.T1.9"  ,  
                    "Ex45.T1.3" ,"Ex45.T1.6", "Ex45.T1.9" ,
                    "HP.T1.3",   "HP.T1.6"  ,"HP.T1.9" ,
                    "RC.T1.3" ,  "RC.T1.6" ,  "RC.T1.9" ,
                    "RP.T1.3" ,  "RP.T1.6" ,  "RP.T1.9"  ,
                    "SC.T1.3" ,  "SC.T1.6"  , "SC.T1.9"  ,
                    "UAB.T1.3" , "UAB.T1.6" ,"UAB.T1.9" )
colnames(kegg_counts_rich) <- names_columnas 


# calculate richness

md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")

summary(colSums(kegg_counts_rich[2:25]))

md$p_func_rich <- specnumber(rrarefy(t(kegg_counts_rich[2:25]), sample = 55297))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "p_func_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(p_func_rich),
                   p_func_rich = mean(p_func_rich))


ggplot(md, aes(x = urban.natural, y = p_func_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  p_func_rich-sd, ymax =  p_func_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("Plasmid functional richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_functional_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(p_func_rich ~ urban.natural + (1|site), data = md), type = 3)
model_log <- lmer(log(p_func_rich) ~ urban.natural + (1|site), data = md)
Anova(lmer(log(p_func_rich) ~ urban.natural + (1|site), data = md))
shapiro.test(residuals(model_log))


#######################
# RPKM
#######################

# 1. prepare the length and the total reads to left join
# import and format the table of lengths per gene
plasmids_len <- all_summary_filtered %>%
  select(c("seq_name", "length"))

colnames(plasmids_len) <- c("Contig", "len")

# select the total reads per sample from the metadata
total_reads <- md %>%
  dplyr::select(sample, total.reads)

# 2. pivot longer the count table to make every gene in every sample a row

kegg_counts_long<- kegg_counts %>%
  ungroup() %>%
  select(-c("kegg_id" , "target_name" ,"target_accession" , "query_name" ,"query_accession","e_value","score" , "bias" , "best_domain_e_value" ,
            "best_domain_score" ,"best_domain_bias", "...1" , "level1" ,"level3" , "KO_description", "seq_name", "conjugation_genes",
            "amr_genes")) %>%
  dplyr::distinct(KO, Contig, .keep_all = TRUE) %>%
  #select(-KO) %>%
  rownames_to_column(var = "row_id") %>% # this is to avoid issues from Contigs and KO_descriptions being the same across rows
  pivot_longer(cols = -c("Contig", "level2", "row_id", "KO"), values_to = "Count", names_to= "sample") %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = ""))

# 3. join the total reads to each row based on sample
kegg_counts_long_plusreads <- left_join(kegg_counts_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
kegg_counts_long_plusreads_plustlen <- left_join(kegg_counts_long_plusreads, plasmids_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_kegg_final <- kegg_counts_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Count", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) 

# make sure all columns are the right class
glimpse(rpkm_kegg_final)


# Sum the counts of each function to run maaslin2
rpkm_kodescription_final <- rpkm_kegg_final %>%
  select(-c("Contig", "row_id", "KO")) %>%
  pivot_longer(-level2, names_to = "sample", values_to = "counts") %>%
  group_by(level2, sample) %>%
  summarize(counts = sum(counts)) %>%
  pivot_wider(names_from = "sample", values_from = "counts") %>%
  drop_na() %>%
  column_to_rownames(var = "level2")


############################
# beta-diversity
############################
rpkm_kegg_final_for_nmds <- rpkm_kegg_final %>%
  select(-c("row_id", "level2", "Contig")) %>%
  distinct(KO, .keep_all = T) %>%
  column_to_rownames(var = "KO")

kegg.bray <- vegdist(t(rpkm_kodescription_final), method="bray")
kegg.bray <- vegdist(t(rpkm_kegg_final_for_nmds), method="bray") # use this one


# Calculate beta dispersion
ptus.disp <- betadisper(kegg.bray, md$urban.natural)
# Test for significant differences in dispersion
permutest(ptus.disp, permutations = 999)
# F N.Perm Pr(>F)
# 0.4343    999  0.531

# Extract distances to centroid for each sample
ptus.disp$distances
# Get group distances to centroid (mean dispersion per group)
ptus.disp$group.distances
# Visualize dispersion
library(patchwork)  # or use gridExtra

# Extract data for ggplot
disp_df <- data.frame(
  Distance = ptus.disp$distances,
  Group = ptus.disp$group
)

p1 <-  ggplot(disp_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = pal2) +
  theme_bw() +
  labs(title = "Distance to Centroid", y = "Distance to centroid") +
  theme(legend.position = "none")

# For the PCoA plot, extract coordinates from the betadisper object
# The vectors are stored in ptus.disp$vectors
pcoa_df <- data.frame(
  PCoA1 = ptus.disp$vectors[, 1],
  PCoA2 = ptus.disp$vectors[, 2],
  Group = ptus.disp$group
)

# Get centroids for each group
centroids_df <- data.frame(
  PCoA1 = ptus.disp$centroids[, 1],
  PCoA2 = ptus.disp$centroids[, 2],
  Group = rownames(ptus.disp$centroids)
)

p2 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse() +
  geom_point(data = centroids_df, aes(x = PCoA1, y = PCoA2), 
             size = 5, shape = 17) +  # Add centroids as triangles
  scale_color_manual(values = pal2) +
  theme_bw() +
  labs(title = "PCoA of Dispersion")

p2 + p1 + plot_layout(widths = c(1.2, 1))

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/kegg_comp_dispersion.pdf", 
       width = 10, height = 4, units = "in")

#ordination (non-multidimensional scaling)
mges.nmds <- metaMDS(kegg.bray, k=2, try = 100)
md$Axis01 = mges.nmds$points[,1]
md$Axis02 = mges.nmds$points[,2]
mges.nmds$stress #0.03989818 the smaller the better, good <0.3)

ggplot(md, aes(Axis01, Axis02))+
  geom_point(aes(color=urban.natural), size=4)+
  stat_ellipse(aes(color = urban.natural))+
  scale_color_manual(values=pal2)+
  theme_bw()+
  theme(legend.position="botoom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_kegg_profile.pdf", device = "pdf", width = 5, height = 5 , units = "in")

adonis2(kegg.bray  ~ urban.natural, data = md, permutations = 999, method = "bray", strata = md$site)
adonis2(kegg.bray  ~ urban.natural, data =  md, permutations = 999, method = "bray") 

#Differential abundance of functions (Maaslin2)
# features as columns and samples as rows
# rownames needs to match the metadata

md_maaslin <- md %>%
  column_to_rownames(var = "sample")

library(Maaslin2)
fit_data = Maaslin2(
  input_data = rpkm_kodescription_final, 
  input_metadata = md_maaslin, 
  output = "/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/", 
  fixed_effects = "urban.natural", 
  transform = "LOG")


sig_families <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/significant_results.tsv", sep = "\t", header = T)
# filter and order

sig_families <- sig_families %>%
  filter(pval < 0.05) %>%
  mutate(color_code = case_when(coef > 0 ~ "Green", 
                           TRUE ~ "Orange")) %>%
  arrange(desc(abs(coef))) %>%
  slice_head(n = 50)

ggplot(sig_families, aes(x = reorder(feature, -coef), y = coef)) + 
  geom_col(aes(fill = color_code)) + 
  scale_fill_manual(values =c("#77966D","#56282D"))+
  theme_bw()+
  xlab("Function")+
  ylab("")+
  coord_flip()+
  theme(legend.position = "none")
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/kegg_level2_maaslin_plasmids.pdf", device = "pdf", width = 11, height = 7 , units = "in")


############################
# abundance of plasmids with antibiotic resistance
############################

rpkm_amr <- rpkm_kegg_final %>%
  filter(level2 %in% c("09111 Xenobiotics biodegradation and metabolism")) %>%
  distinct(level2, Contig, .keep_all = TRUE) %>%
  select(-c("row_id", "level2", "KO")) %>%
  column_to_rownames("Contig")

colSums(rpkm_amr)

md$xenobio_relab <- colSums(rpkm_amr)
  
# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "xenobio_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(xenobio_relab),
                   xenobio_relab= mean(xenobio_relab))

pal1 = c("#FF823B", "#359951")

p <- ggplot(md, aes(x = urban.natural, y = xenobio_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  xenobio_relab-sd, ymax =  xenobio_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("Xenobiotic resistance relative abundance")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

p_marg

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/xenobiotic_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")

#ggsave("/Volumes/BunnyBike/mge_urban/local/figures/xenobiotic_resistance_relab.pdf", device = "pdf", width = 5, height = 5 , units = "in")
model <- lmer(xenobio_relab ~ urban.natural + (1|site), data = md)
Anova(model, type = 3)
shapiro.test(residuals(model))
model_log <- lmer(log(xenobio_relab) ~ urban.natural + (1|site), data = md)
Anova(lmer(log(xenobio_relab) ~ urban.natural + (1|site), data = md))
shapiro.test(residuals(model_log))

# richness

kegg_counts_amr <- kegg_counts %>%
  filter(level2 %in% c("09111 Xenobiotics biodegradation and metabolism")) %>%
  ungroup() %>%
  select(-c("kegg_id" , "target_name" ,"target_accession" , "query_name" ,"query_accession","e_value","score" , "bias" , "best_domain_e_value" ,
            "best_domain_score" ,"best_domain_bias", "...1" , "level1" ,"level3" , "level2", "KO", "KO_description", "seq_name", "conjugation_genes",
            "amr_genes")) %>%
  distinct(Contig, .keep_all = T) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  pivot_wider(names_from = "sample", values_from= "counts") %>%
  column_to_rownames(var = "Contig")
  

summary(colSums(kegg_counts_amr))

md$xen_rich  <- specnumber(rrarefy(t(kegg_counts_amr), sample = 1305))
# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "xen_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(xen_rich),
                   xen_rich = mean(xen_rich))

pal1 = c("#FF823B", "#359951")

ggplot(md, aes(x = urban.natural, y = xen_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  xen_rich-sd, ymax = xen_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("Number of plasmids with Xenobiotic resistance genes")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/kegg_xenres_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")





# First let's prepare the rpkm table here so we can easily modify it for each group
#######################
# RPKM
#######################
# Level 1: 
# we don't need all the steps, just these: 
kegg_counts_long<- kegg_counts %>%
  ungroup() %>%
  select(-c("kegg_id" , "target_name" ,"target_accession" , "query_name" ,"query_accession","e_value","score" , "bias" , "best_domain_e_value" ,
            "best_domain_score" ,"best_domain_bias", "...1" , "level2" , "level3" ,  "KO_description", "seq_name", "conjugation_genes",
            "amr_genes")) %>%
  dplyr::distinct(KO, Contig, .keep_all = TRUE) %>%
  #select(-KO) %>%
  rownames_to_column(var = "row_id") %>% # this is to avoid issues from Contigs and KO_descriptions being the same across rows
  pivot_longer(cols = -c("Contig", "level1", "row_id", "KO"), values_to = "Count", names_to= "sample") %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = ""))

# 3. join the total reads to each row based on sample
kegg_counts_long_plusreads <- left_join(kegg_counts_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
kegg_counts_long_plusreads_plustlen <- left_join(kegg_counts_long_plusreads, plasmids_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_kegg_final_level1 <- kegg_counts_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Count", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) 

# level 2:
# that is the original rpkm_counts_final
rpkm_kegg_final_level2 <- rpkm_kegg_final

# Level 3: 
kegg_counts_long<- kegg_counts %>%
  ungroup() %>%
  select(-c("kegg_id" , "target_name" ,"target_accession" , "query_name" ,"query_accession","e_value","score" , "bias" , "best_domain_e_value" ,
            "best_domain_score" ,"best_domain_bias", "...1" , "level2" , "level1" ,  "KO_description", "seq_name", "conjugation_genes",
            "amr_genes")) %>%
  dplyr::distinct(KO, Contig, .keep_all = TRUE) %>%
  #select(-KO) %>%
  rownames_to_column(var = "row_id") %>% # this is to avoid issues from Contigs and KO_descriptions being the same across rows
  pivot_longer(cols = -c("Contig", "level3", "row_id", "KO"), values_to = "Count", names_to= "sample") %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = ""))

# 3. join the total reads to each row based on sample
kegg_counts_long_plusreads <- left_join(kegg_counts_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
kegg_counts_long_plusreads_plustlen <- left_join(kegg_counts_long_plusreads, plasmids_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_kegg_final_level3 <- kegg_counts_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Count", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) 



####### Motility
kegg_counts_group <- rpkm_kegg_final_level1 %>%
  filter(KO %in% c("K22221","K02282", "K02282","K02424")) %>% # motility
  distinct(level1, Contig, .keep_all = TRUE) %>%
  select(-c("row_id", "level1", "KO")) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(counts = sum(counts))

####### Toxins
kegg_counts_group<- rpkm_kegg_final_level3 %>%
  filter(str_detect(level3, "toxin")) %>% # toxins
  distinct(level3, Contig, .keep_all = TRUE) %>%
  select(-c("row_id", "level3", "KO")) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(counts = sum(counts))

####### Defense
kegg_counts_group <- rpkm_kegg_final_level3 %>%
  filter(str_detect(level3, "defense")) %>%
  distinct(level3, Contig, .keep_all = TRUE) %>%
  select(-c("row_id", "level3", "KO")) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(counts = sum(counts))

####### Metabolism
kegg_counts_group <- rpkm_kegg_final_level1 %>%
  filter(str_detect(level1, "Metabolism")) %>%
  distinct(level1, Contig, .keep_all = TRUE) %>%
  select(-c("row_id", "level1", "KO")) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(counts = sum(counts))

summary(colSums(kegg_counts_group))

md_mot <- md %>%
  left_join(., kegg_counts_group, by = "sample")

# Line with error bar
mrelab_lineplot <- md_mot %>%
  dplyr::select(c("urban.natural", "counts")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(counts),
                   counts = mean(counts))

p <- ggplot(md_mot, aes(x = urban.natural, y = counts)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  geom_pointrange(aes(ymin =  counts-sd, ymax = counts+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("Motility gene abundance (RPKM)")+
  ylab("Bacterial toxins gene abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

p_marg

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/motility_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")


Anova(lmer(counts ~ urban.natural + (1|site), data = md_mot), type = 3)
model <- lmer(counts ~ urban.natural + (1|site), data = md_mot)
shapiro.test(residuals(model))

model_log <- lmer(log(counts) ~ urban.natural + (1|site), data = md_mot)
Anova(lmer(log(counts) ~ urban.natural + (1|site), data = md_mot))
shapiro.test(residuals(model_log))

# Motility
#Chisq Df Pr(>Chisq)  
#(Intercept)   0.0125  1    0.91115  
#urban.natural 4.3082  1    0.03793 *
#log: urban.natural 4.0637  1    0.04381 *

# Toxins
#Chisq Df Pr(>Chisq)  
#(Intercept)   3.5288  1    0.06031 
#urban.natural 2.0519  1    0.15201 

# Defense
#Chisq Df Pr(>Chisq)  
#(Intercept)   1.8413  1    0.17479  
#urban.natural 5.5366  1    0.01862 *
# log:urban.natural 8.5866  1   0.003386

# Metabolism
#Chisq Df Pr(>Chisq)    
#(Intercept)   24.3693  1  7.952e-07 ***
#  urban.natural  2.9351  1    0.08668 .


############################################################### ptu counts









# MODEL DIAGNOSTICS

# models
# p_gene_rich
model <- lmer(p_gene_rich ~ urban.natural + (1|site), data = md)
# p_func_rich
model <- lmer(p_func_rich ~ urban.natural + (1|site), data = md)
# xenobio_relab
model <- lmer(xenobio_relab ~ urban.natural + (1|site), data = md)
# xenobio_relab
model <- lmer(motility_counts ~ urban.natural + (1|site), data = md_mot)

## Model diagnostics
# Residual plots
library(patchwork)

# Extract model data
model_df <- data.frame(
  Fitted = fitted(model),
  Residuals = residuals(model)
)

# Residual plot
p1 <- ggplot(model_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(title = "Residuals vs Fitted", 
       x = "Fitted values", 
       y = "Residuals")

# Q-Q plot
p2 <- ggplot(model_df, aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_bw() +
  labs(title = "Normal Q-Q Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Combine and save
p1 + p2

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/p_func_rich_model_diagnostics.pdf", 
       width = 10, height = 4, units = "in")

# Shapiro-Wilk test for normality
shapiro.test(residuals(model))



# RESULTS: 

# p_gene_rich : W = 0.94741, p-value = 0.2379
# normality not violated

# p_func_rich: W = 0.91628, p-value = 0.04835     !!
model_log <- lmer(log(p_func_rich) ~ urban.natural + (1|site), data = md)
Anova(lmer(log(p_func_rich) ~ urban.natural + (1|site), data = md))
shapiro.test(residuals(model_log))
# Square root transformation
model_sqrt <- lmer(sqrt(p_func_rich) ~ urban.natural + (1|site), data = md)
Anova(lmer(sqrt(p_func_rich) ~ urban.natural + (1|site), data = md))
shapiro.test(residuals(model_sqrt))
glmer(p_func_rich ~ urban.natural + (1|site), 
      data = md, family = poisson)

# Permutation test (doesn't assume normality)
model_perm <- lmp(p_func_rich ~ urban.natural, data = md)
summary(model_perm)
#The QQ plot looks ok so we can proceed with caution.
#Both log transformation and sqrt transformation made normality worst, 
# but all of them give me significant results. Permutation test also gives significant results. 




check_model(model, check=c("linearity", "homogeneity", "normality"))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/p_gene_rich_model_diagnostics_performace.pdf", 
       width = 6, height = 4, units = "in")
