# analysis of bacmet results for urban paper

library(tidyverse)


ptus_bacmet <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_bacmet.out", sep = "\t", header = F)

colnames(ptus_bacmet) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                                "qend", "sstart", "send", "evalue", "bitscore")
ptus_bacmet <- ptus_bacmet %>%
  distinct(qseqid, .keep_all = T)

ptus_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_count_table.txt", sep = "\t", header = T)

md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")


#####################
# Counts
####################

count_bacmet_ptus <-ptus_counts %>%
  filter(Contig %in% ptus_bacmet$qseqid) %>%
  column_to_rownames(var = "Contig")

summary(colSums(count_bacmet_ptus ))

md$bacmet_rich  <- specnumber(rrarefy(t(count_bacmet_ptus), sample = 4087))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "bacmet_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(bacmet_rich),
                   bacmet_rich = mean(bacmet_rich))

pal1 = c("#FF823B", "#359951")
pal1 = c("#56282D", "#77966D")

ggplot(md, aes(x = urban.natural, y = bacmet_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  bacmet_rich-sd, ymax = bacmet_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("HMRG plasmid richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/bacmet_richness.pdf", device = "pdf", width = 5, height = 5 , units = "in")




#############
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

count_bacmet_plasmids_long <- count_bacmet_ptus %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "Counts") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) 

# 3. join the total reads to each row based on sample
bacmet_long_plusreads <- left_join(count_bacmet_plasmids_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
bacmet_long_plusreads_plustlen <- left_join(bacmet_long_plusreads, ptus_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_bacmet_final <- bacmet_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Counts*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Counts", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig") %>%
  drop_na()


############################
# abundance
############################
sum_rel_ab_bacmet <- as.data.frame(t(rpkm_bacmet_final)) 

sum_rel_ab_bacmet_final <- sum_rel_ab_bacmet %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(bacmet_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(bacmet_relab) %>%
  rownames_to_column(var = "sample")

md_final<- md %>%
  left_join(., sum_rel_ab_bacmet_final, by = "sample")


# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "bacmet_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(bacmet_relab),
                   bacmet_relab  = mean(bacmet_relab))

p <- ggplot(md_final, aes(x = urban.natural, y = bacmet_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =bacmet_relab-sd, ymax = bacmet_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("HMRG abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))

# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

p_marg

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/hmrgs_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(bacmet_relab ~ urban.natural + (1|site), data = md_final), type = 3)
model <- lmer(bacmet_relab ~ urban.natural + (1|site), data = md_final)
shapiro.test(residuals(model))

model_log <- lmer(log(bacmet_relab) ~ urban.natural + (1|site), data = md_final)
Anova(lmer(log(bacmet_relab) ~ urban.natural + (1|site), data = md_final))
shapiro.test(residuals(model_log))


#ggsave
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/bacmetfinder_relab_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")



######################### Types of metal resistance genes
# categorize
mapping_file <- read.table("/Volumes/BunnyBike/landuse_functionaltraits/metals/BacMet2_EXP.753.mapping.txt", sep = "\t", header = T)

mapping_file_join <- mapping_file %>%
  select(c("BacMet_ID", "Compound"))

colnames(mapping_file)
# join mapping to ptu counts
ptus_bacmet_metals <- ptus_bacmet %>%
  separate(sseqid, sep = "\\|", into=c("BacMet_ID", "Gene_name", "sm1", "sm2", "sm3")) %>%
  select(-c("sm1", "sm2", "sm3")) %>%
  left_join(mapping_file_join, by = "BacMet_ID")



# plot the number of genes of each type overall
ptus_bacmet_metals_count <- ptus_bacmet_metals %>%
  group_by(Compound) %>%
  summarize(n = n()) %>%
  separate(Compound, sep = ", ", into = c("C1", "C2", "C3", "C4", "C5")) %>%
  pivot_longer(
    cols = starts_with("C"),
    names_to = "column",
    values_to = "metal") %>%
  filter(!is.na(metal)) %>%
  select(-column) %>%
  mutate(metal = gsub(metal, pattern = " \\[class:.*\\]$", replacement = "")) %>%
  group_by(metal) %>%
  summarize(n = sum(n))
  

ptus_bacmet_metals_count$metal
ggplot(ptus_bacmet_metals_count, aes (x = reorder(metal, n), y = n))+
  geom_col()+
  coord_flip()+
  xlab("Number of genes")+
  ylab("")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/metals_number_cols_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")



# plot the number of genes in urban vs natural of each type
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

ptus_bacmet_metals_join <- ptus_bacmet_metals %>%
  select(c("qseqid", "Compound"))

colnames(ptus_bacmet_metals_join) <- c("plasmid", "class")  

ptus_bacmet_metals_join_class <- ptus_bacmet_metals_join %>%
  left_join(., ptu_pab, by = "plasmid") %>%
  drop_na() %>%
  select(-plasmid) %>%
  pivot_longer(-class, names_to = "sample", values_to="n") %>%
  separate(class, sep = ", ", into = c("C1", "C2", "C3", "C4", "C5")) %>%
  pivot_longer(
    cols = starts_with("C"),
    names_to = "column",
    values_to = "metal") %>%
  filter(!is.na(metal)) %>%
  select(-column) %>%
  group_by(metal, sample) %>%
  summarize(counts = sum(n)) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural")) %>%
  filter(metal %in% c("Copper (CU)", "Nickel (Ni)", "Iron (Fe)", "Tungsten (W)", 
                      "Cobalt (Co)", "Arsenic (As)", "Zinc (Zn)", "Mercury (Hg)",
                      "Lead (Pb)")) 
  
  ptus_bacmet_metals_join_class$metal

ggplot(ptus_bacmet_metals_join_class, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~metal)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number of genes")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/metals_number_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")
