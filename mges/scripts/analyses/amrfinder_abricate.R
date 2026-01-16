# analysis of results from AMRfinder urban plasmids project

library(tidyverse)
library(vegan)
library(ggExtra)


############################################################################# AMR FINDER
md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")

###########
##### AMR in plasmids
###########

# load summary plasmids (info on their AMR finder here)
all_summary_filtered <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/circular_10kb_plasmids_summary.csv")


amr_plasmids <- all_summary_filtered %>%
  filter(!is.na(amr_genes))
# 26/213

amr_plasmids_conj <- all_summary_filtered %>%
  filter(!is.na(amr_genes)) %>%
  filter(!is.na(conjugation_genes))
# 1

conj_plasmids <- all_summary_filtered %>%
  filter(!is.na(conjugation_genes))

###########
##### AMR in all contigs
###########
amr_all <- read.table("/Volumes/BunnyBike/mge_urban/local/amrfinder/all_contigs_amrfinder.out", sep = "\t", header = T)

# remove the header lines
amr_all <- amr_all %>%
  filter(!Start == "Start") #%>%
  #distinct(Contig.id)
# 643



#####################
# Counts
####################
plasmid_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_count_table.txt", sep = "\t", header = T)

count_amr_plasmids <- plasmid_counts %>%
  filter(Contig %in% amr_plasmids$seq_name) %>%
  column_to_rownames(var = "Contig")

summary(colSums(count_amr_plasmids))

md$amr_rich  <- specnumber(rrarefy(t(count_amr_plasmids), sample = 633))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "amr_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(amr_rich),
                   amr_rich = mean(amr_rich))

pal1 = c("#FF823B", "#359951")

ggplot(md, aes(x = urban.natural, y = amr_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  amr_rich-sd, ymax = amr_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("AMR plasmid richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/amrfinder_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")



##############
# RPKM
##############

all_summary_filtered <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/circular_10kb_plasmids_summary.csv")
# 1. prepare the length and the total reads to left join
# import and format the table of lengths per gene
ptus_len <- all_summary_filtered %>%
  select(c("seq_name", "length"))

colnames(ptus_len) <- c("Contig", "len")

# select the total reads per sample from the metadata
total_reads <- md %>%
  dplyr::select(sample, total.reads)

count_amr_plasmids_long <- count_amr_plasmids %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "Counts") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) 

# 3. join the total reads to each row based on sample
amr_long_plusreads <- left_join(count_amr_plasmids_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
amr_long_plusreads_plustlen <- left_join(amr_long_plusreads, ptus_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_amr_final <- amr_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Counts*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Counts", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig") %>%
  drop_na()

############################
# abundance
############################
sum_rel_ab_amr <- as.data.frame(t(rpkm_amr_final)) 

sum_rel_ab_amr_final <- sum_rel_ab_amr %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(amr_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(amr_relab) %>%
  rownames_to_column(var = "sample")

md_final<- md %>%
  left_join(., sum_rel_ab_amr_final, by = "sample")


# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "amr_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(amr_relab),
                   amr_relab  = mean(amr_relab))

p <- ggplot(md_final, aes(x = urban.natural, y = amr_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =amr_relab-sd, ymax = amr_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("AMR abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))

# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

p_marg

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/AMR_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")

#ggsave
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/amrfinder_relab_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(amr_relab ~ urban.natural + (1|site), data = md_final), type = 3)
model <- lmer(amr_relab ~ urban.natural + (1|site), data = md_final)
shapiro.test(residuals(model))

model_log <- lmer(log(amr_relab) ~ urban.natural + (1|site), data = md_final)
Anova(lmer(log(amr_relab) ~ urban.natural + (1|site), data = md_final))
shapiro.test(residuals(model_log))




############################################################################# ABRICATE
###########
##### Abricate all plasmids
###########
plasmids_abricate <- read_table("/Volumes/BunnyBike/mge_urban/local/abricate/plasmids/ptus_abricate.tab")

mges_abricate <- read_table("/Volumes/BunnyBike/mge_urban/local/abricate/plasmids/mgtus_abricate.tab")

# 
all_abricate <- read_table("/Volumes/BunnyBike/mge_urban/local/abricate/all_abricate.tab")


all_abricate <- all_abricate %>%
  filter(!SEQUENCE == "SEQUENCE")






########################################################################## CONJUGATIVE PLASMIDS
count_conj_plasmids <- plasmid_counts %>%
  filter(Contig %in% conj_plasmids$seq_name) %>%
  column_to_rownames(var = "Contig")


summary(colSums(count_conj_plasmids))

md$conj_rich  <- specnumber(rrarefy(t(count_conj_plasmids), sample = 15))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "conj_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(conj_rich),
                   conj_rich = mean(conj_rich))

pal1 = c("#FF823B", "#359951")

ggplot(md, aes(x = urban.natural, y = conj_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  conj_rich-sd, ymax = conj_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Conjugative plasmid richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/conjugative_plasmid_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")


##############
# RPKM
##############

all_summary_filtered <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/circular_10kb_plasmids_summary.csv")
# 1. prepare the length and the total reads to left join
# import and format the table of lengths per gene
ptus_len <- all_summary_filtered %>%
  select(c("seq_name", "length"))

colnames(ptus_len) <- c("Contig", "len")

# select the total reads per sample from the metadata
total_reads <- md %>%
  dplyr::select(sample, total.reads)

count_conj_plasmids_long <- count_conj_plasmids %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "Counts") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) 

# 3. join the total reads to each row based on sample
conj_long_plusreads <- left_join(count_conj_plasmids_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
conj_long_plusreads_plustlen <- left_join(conj_long_plusreads, ptus_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_conj_final <- conj_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Counts*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Counts", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig") %>%
  drop_na()

############################
# abundance
############################
sum_rel_ab_conj <- as.data.frame(t(rpkm_conj_final)) 

sum_rel_ab_conj_final <- sum_rel_ab_conj %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(conj_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(conj_relab) %>%
  rownames_to_column(var = "sample")

md_final<- md %>%
  left_join(., sum_rel_ab_conj_final, by = "sample")


# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "conj_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(conj_relab),
                   conj_relab  = mean(conj_relab))

ggplot(md_final, aes(x = urban.natural, y = conj_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =conj_relab-sd, ymax = conj_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Conjugative plasmid abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/conjplasmid_relab_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")




## Types of AMR, from AMRFinder

amr_db <- read_tsv("/Volumes/BunnyBike/mge_urban/local/amrfinder/NCBIfam-AMRFinder.tsv")

amr_db <- amr_db %>%
  mutate(amr_genes = gsub(hmm_accession, pattern = ".[0-9]*$", replacement = ""))

amr_types <- amr_plasmids  %>%
  select(c("seq_name", "amr_genes")) %>%
  separate_rows(amr_genes, sep = ";") %>%
  left_join(., amr_db, by = "amr_genes")

amr_class <- amr_types %>%
  group_by(gene_symbol) %>%
  summarize(count = n())

ggplot(amr_class, aes(x = reorder(gene_symbol, count), y = count))+
  geom_col() +
  xlab("")+
  ylab("Number of genes")+
  theme_bw()+
  coord_flip()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/amrfinder_gene_symbol.pdf", device = "pdf", width = 5, height = 5 , units = "in")




# Types amrs presence in urban vs. natural
# join this information with presence absence in urban and natural

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

amr_counts <- amr_types %>%
  select(c("seq_name", "class"))

colnames(amr_counts) <- c("plasmid", "class")
  
amr_counts_class <-amr_counts  %>%
  left_join(., ptu_pab, by = "plasmid") %>%
  drop_na() %>%
  select(-plasmid) %>%
  pivot_longer(-class, names_to = "sample", values_to="count") %>%
  group_by(class, sample) %>%
  summarize(counts = sum(count)) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                                               "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                                                 TRUE ~ "natural"))

ggplot(amr_counts_class, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~class)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number of genes")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/amr_class_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")



# Check abundance, not presence of genes

# NOTE: RUN THE RPKM CALCULATION FROM ptus_analyses.R
ptu_counts_join <- rpkm_ptus_final %>%
  rownames_to_column(var = "plasmid") %>%
  pivot_longer(-c("plasmid"), names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "-", replacement = "\\.")) %>%
  pivot_wider(names_from = "sample", values_from = "count")
  
  

amr_counts_class <- amr_counts  %>%
  left_join(., ptu_counts_join, by = "plasmid") %>%
  drop_na() %>%
  select(-plasmid) %>%
  pivot_longer(-class, names_to = "sample", values_to="count") %>%
  group_by(class, sample) %>%
  summarize(counts = sum(count)) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural"))

ggplot(amr_counts_class, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~class)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/amr_class_rpkm_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")


  