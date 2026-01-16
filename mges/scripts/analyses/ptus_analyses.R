# analysis of PTU count table for URBAN plasmids project


library(tidyverse)
library(vegan)
library(see)
library(lme4)
library(car)
library(performance)
library(lmPerm)
library(patchwork)
pal1 = c("#56282D", "#77966D")


# load count table
ptu_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_count_table.txt", sep = "\t", header = T)


# load metadata 
md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")

md <- md %>%
  mutate(sample = gsub(sample, pattern = "\\.", replacement = "-"))


# Match PTU and md tables
ptu_counts_long <- ptu_counts %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "Counts") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "\\.", replacement = "-")) 


ptu_counts_filtered <- ptu_counts_long %>%
  filter(sample %in% md$sample) %>%
  pivot_wider(names_from = Contig, values_from = Counts ) %>%
  arrange(desc(sample)) %>%
  column_to_rownames(var = "sample") 


md_final<- md %>%
  arrange(desc(sample)) %>%   # check that the
  column_to_rownames(var = "sample") 

rownames(ptu_counts_filtered) == rownames(md_final)



# Richness
summary(rowSums(ptu_counts_filtered))

md_final$ptu_rich <- specnumber(rrarefy(ptu_counts_filtered, sample = 5813))
md_final$ptu_shannon <- diversity(rrarefy(ptu_counts_filtered, sample = 5813))
md_final$ptu_evenness <- md_final$ptu_shannon/log(md_final$ptu_rich)

# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "ptu_evenness")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(ptu_evenness),
                   ptu_evenness  = mean(ptu_evenness))

ggplot(md_final, aes(x = urban.natural, y = ptu_evenness)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin = ptu_evenness-sd, ymax = ptu_evenness+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("PTU Evenness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptu_evenness_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# STATS
model <- lmer(ptu_evenness ~ urban.natural + (1|site), data = md_final)
model <- lmer(ptu_rich ~ urban.natural + (1|site), data = md_final)
model <- lmer(ptu_shannon ~ urban.natural + (1|site), data = md_final)


model = lmer(ptus_relab ~ urban.natural + (1|site), data = md_final)

Anova(model, type = 3)
#Response: ptu_rich
#Chisq Df Pr(>Chisq)    
#(Intercept)   1594.828  1  < 2.2e-16 ***
#  urban.natural   16.244  1  5.569e-05 ***

#Response: ptu_shannon
#Chisq Df Pr(>Chisq)    
#(Intercept)   1118.3296  1     <2e-16 ***
#  urban.natural    5.1538  1     0.0232 *

#Response: ptu_evenness
#Chisq Df Pr(>Chisq)    
#(Intercept)   1375.2326  1    < 2e-16 ***
#  urban.natural    3.2635  1    0.07084 .


# lm? 
Anova(lm(ptu_rich ~ urban.natural * vegetation_structure, data = md_final))


## Model diagnostics
# Residual plots

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

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptu_relab_model_diagnostics.pdf", 
       width = 10, height = 4, units = "in")

# Shapiro-Wilk test for normality
shapiro.test(residuals(model))

# Check model (performance)
check_model(model, check=c("linearity", "homogeneity", "normality"))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptu_relab_checkmodel.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# Evenness and Shannon are not normal: 
# and they are continuous data, so no glmm with poison
# Log transform
model_log <- lmer(log(ptus_relab) ~ urban.natural + (1|site), data = md_final)
Anova(lmer(log(ptus_relab) ~ urban.natural + (1|site), data = md_final))
shapiro.test(residuals(model_log))
# shannon
# anova: urban.natural 5.0008  1    0.02534 *
# shapiro: W = 0.86515, p-value = 0.004234
# relab
# anova: urban.natural 13.497  1   0.000239 ***
# shapiro: W = 0.90966, p-value = 0.03468 
# log transform doesn't help

# Square root transformation
model_sqrt <- lmer(sqrt(ptus_relab) ~ urban.natural + (1|site), data = md_final)
Anova(lmer(sqrt(ptus_relab) ~ urban.natural + (1|site), data = md_final))
shapiro.test(residuals(model_sqrt))
# W = 0.8556, p-value = 0.002783
# nor does sqrt

# GLMM (do not use for relative abundance): 
shannon.glmm <- glmer(ptu_shannon ~ urban.natural + (1|site), data = md_final, family = Gamma(link = "log"))
# Identity link
model_identity <- glmer(ptu_shannon ~ urban.natural + (1|site), 
                        data = md_final, 
                        family = Gamma(link = "identity"))

# Inverse link (default for Gamma)
model_inverse <- glmer(ptu_shannon ~ urban.natural + (1|site), 
                       data = md_final, 
                       family = Gamma(link = "inverse"))

library(DHARMa)
simulateResiduals(model_inverse, plot = T)
# permutation test: 
model_perm <- lmp(ptus_relab ~ urban.natural, data = md_final)
summary(model_perm)


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

# 3. join the total reads to each row based on sample
ptus_long_plusreads <- left_join(ptu_counts_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
ptus_long_plusreads_plustlen <- left_join(ptus_long_plusreads, ptus_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_ptus_final <- ptus_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Counts*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Counts", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig") %>%
  drop_na()



############################
# abundance
############################
sum_rel_ab_ptus <- as.data.frame(t(rpkm_ptus_final)) 

sum_rel_ab_ptus_final <- sum_rel_ab_ptus %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(ptus_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(ptus_relab) %>%
  rownames_to_column(var = "sample")

md_final <- md_final %>%
  rownames_to_column(var = "sample") %>%
  left_join(., sum_rel_ab_ptus_final, by = "sample")


# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "ptus_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(ptus_relab),
                   ptus_relab  = mean(ptus_relab))

ggplot(md_final, aes(x = urban.natural, y = ptus_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin = ptus_relab-sd, ymax = ptus_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("PTUS abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_relab_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")

model = lmer(ptus_relab ~ urban.natural + (1|site), data = md_final)
Anova(model, type = 3)
Anova(lm(ptus_relab ~ urban.natural * vegetation_structure, data = md_final))


## Model diagnostics above


############################
# beta-diversity
############################
ptus.bray <- vegdist(t(rpkm_ptus_final), method="bray")

# Calculate beta dispersion
ptus.disp <- betadisper(ptus.bray, md_final$urban.natural)
# Test for significant differences in dispersion
permutest(ptus.disp, permutations = 999)
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
  scale_fill_manual(values = pal1) +
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
  scale_color_manual(values = pal1) +
  theme_bw() +
  labs(title = "PCoA of Dispersion")

p2 + p1 + plot_layout(widths = c(1.2, 1))

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_comp_dispersion.pdf", 
       width = 10, height = 4, units = "in")


# Anova for testing
anova(ptus.disp)



#ordination (non-multidimensional scaling)
ptus.nmds <- metaMDS(ptus.bray, k=2, try = 100)
md_final$Axis01 = ptus.nmds$points[,1]
md_final$Axis02 = ptus.nmds$points[,2]
ptus.nmds$stress #0.03989818 the smaller the better, good <0.3)

# urban.natural
ggplot(md_final, aes(Axis01, Axis02))+
  geom_point(aes(color=landuse2), size=4)+
  #stat_ellipse(aes(color = landuse2))+
  scale_color_manual(values=pal2)+
  theme_bw()+
  theme(legend.position="bottom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_beta_ecosystem.pdf", device = "pdf", width = 5, height = 5 , units = "in")

adonis2(ptus.bray   ~ urban.natural, data = md_final, permutations = 999, method = "bray", strata = md_final$site)
adonis2(ptus.bray   ~ urban.natural, data =  md_final, permutations = 999, method = "bray") 
adonis2(ptus.bray   ~ urban.natural * vegetation_structure, data =  md_final, permutations = 999, method = "bray") 



# vegetation
ggplot(md_final, aes(Axis01, Axis02))+
  geom_point(aes(shape=vegetation_structure, color = urban.natural), size=4)+
  stat_ellipse(aes(color = urban.natural))+
  scale_color_manual(values=pal1)+
  theme_bw()+
  theme(legend.position="bottom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_beta_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")



############################
# length
############################
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

length_plasmids <- all_summary_filtered %>%
  select(c("seq_name", "length"))

colnames(length_plasmids) <- c("plasmid", "length")



# Length by site

length_sample_for_site <- ptu_pab  %>%
  pivot_longer(-plasmid, names_to = "sample", values_to = "pab") %>%
  left_join(., length_plasmids , by = "plasmid") %>%
  filter(!pab == 0) %>%
  drop_na() 

site_sample <- md %>%
  select(sample, site)

length_sample_for_site_plot <- length_sample_for_site %>%
  mutate(sample = gsub(sample, pattern = "\\.", replacement = "-")) %>%
  left_join(site_sample, by = "sample") %>%
  mutate(urban.natural = case_when(site %in% c("Himmel Park","Reid Park") ~ "Urban",
                                   TRUE ~ "Natural"))

length_sample_for_site_plot$site <- factor(length_sample_for_site_plot$site, levels = c("Himmel Park", "Reid Park", "11" ,  "8" ,  "Ex45", "Rose Canyon" , "Sabino Canyon" ,"UAB" ))

ggplot(length_sample_for_site_plot, aes(x = site, y = length))+
  geom_jitter(aes(color = urban.natural))+
  geom_violin(alpha = 0.4, aes(fill = urban.natural, color = urban.natural)) +
  coord_flip() +
  scale_color_manual(values=pal1)+
  ylab("Plasmid length") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_length_violin.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# barplot: 
md_len_site <- md %>%
  mutate(sample = gsub(sample, pattern = "-", replacement = "\\.")) %>%
  left_join(length_sample_for_site, by = "sample") %>%
  select(c("sample", "urban.natural", "site", "length")) %>%
  group_by(site) %>%
  summarize(mean_length_site = mean(length), 
            sd_site = sd(length)) %>%
  mutate(urban.natural = case_when(site %in% c("Himmel Park","Reid Park") ~ "Urban",
                                   TRUE ~ "Natural"))
  
md_len_site$site <- factor(md_len_site$site, levels = c("Himmel Park", "Reid Park", "11" ,  "8" ,  "Ex45", "Rose Canyon" , "Sabino Canyon" ,"UAB" ))

ggplot(md_len_site, aes(x = site, color = urban.natural, y = mean_length_site)) +
  geom_point()+
  geom_pointrange(aes(ymin =  mean_length_site-sd_site, ymax =  mean_length_site+sd_site, 
                      color = urban.natural),size = 0.5, linewidth = 1, data = md_len_site)+
  geom_line(size = 0.5, alpha = 0.5)+
  coord_flip() +
  scale_color_manual(values=pal1)+
  ylab("Plasmid length") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_length_linerange_site.pdf", device = "pdf", width = 5, height = 5 , units = "in")




# Length boxplot
length_sample <- ptu_pab  %>%
  pivot_longer(-plasmid, names_to = "sample", values_to = "pab") %>%
  left_join(., length_plasmids , by = "plasmid") %>%
  filter(!pab == 0) %>%
  drop_na() %>%
  group_by(sample) %>%
  summarize(mean_length = mean(length), 
            sd_length = sd(length)) %>%
  mutate(urban.natural = case_when(str_detect(sample, pattern = "HP") ~ "Urban",
                                   str_detect(sample, pattern = "RP") ~ "Urban",
                                   TRUE ~ "Natural"))


# linerange per sample
length_sample$sample <- factor(length_sample$sample, levels = c("HP.T1.3","HP.T1.6", "HP.T1.9",
                                                              "RP.T1.3" ,  "RP.T1.6" ,  "RP.T1.9",
                                                              "11B.T1.3" , "11B.T1.6" , "11B.T1.9" ,
                                                              "8.T1.3",    "8.T1.6" ,   "8.T1.9" ,
                                                              "Ex45.T1.3", "Ex45.T1.6" ,"Ex45.T1.9",
                                                              "UAB.T1.3", "UAB.T1.6",  "UAB.T1.9", 
                                                              "RC.T1.3" ,  "RC.T1.6" ,  "RC.T1.9", 
                                                              "SC.T1.3" ,  "SC.T1.6" ,  "SC.T1.9"  ))

ggplot(length_sample, aes(x = sample, color = urban.natural, y = mean_length)) +
  geom_point()+
  geom_pointrange(aes(ymin =  mean_length-sd_length, ymax =  mean_length+sd_length, 
                      color = urban.natural),size = 0.5, linewidth = 1, data = length_sample, alpha = 0.7)+
  #geom_line(size = 0.5, alpha = 0.5)+
  coord_flip() +
  scale_color_manual(values=pal1)+
  ylab("Plasmid length") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_length_linerange_sample.pdf", device = "pdf", width = 5, height = 5 , units = "in")



# Line with error bar
length_sample_tojoin <- length_sample %>%
  select(c("sample", "mean_length")) %>%
  mutate(sample = str_replace(sample, "\\.", "-")) %>%
  mutate(sample = str_replace(sample, "\\.", "-"))

md_length <- md %>%
  left_join(.,length_sample_tojoin, by = "sample")

mrelab_lineplot <- md %>%
  left_join(.,length_sample_tojoin, by = "sample") %>%
  dplyr::select(c("urban.natural", "mean_length")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(mean_length),
                   mean_length = mean(mean_length))

# this doesn't work, fix
ggplot(md, aes(x = urban.natural, y = mean_length)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  mean_length-sd, ymax =  mean_length+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average plasmid length")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_length_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")


Anova(lmer(mean_length ~ urban.natural + (1|site), data = md_length), type = 3)



############################
# Differential abundance
############################

library(Maaslin2)

md_maaslin <- md %>%
  column_to_rownames(var = "sample")

# features as columns and samples as rows
# rownames needs to match the metadata
fit_data = Maaslin2(
  input_data = rpkm_ptus_final, 
  input_metadata = md_maaslin, 
  output = "/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_maaslin", 
  fixed_effects = "urban.natural", 
  random_effects = c("site"),
  transform = "LOG")


sig_families <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_maaslin/significant_results.tsv", sep = "\t", header = T)
# filter and order
sig_families <- sig_families %>%
  filter(qval < 1e-02) %>%
  arrange(desc(abs(coef))) 


# Join conjugation and AMRFinder information
## conjugation
conjugation <- all_summary_filtered %>%
  filter(!is.na(conjugation_genes))

conjugation_tojoin <- conjugation %>%
  select(c("seq_name", "conjugation_genes")) %>%
  mutate(seq_name = gsub(seq_name, pattern = "-", replacement = "\\."))

colnames(conjugation_tojoin) <- c("feature", "conjugation_genes")

# amr finder
amr_finder <- all_summary_filtered %>%
  filter( !is.na((amr_genes)))

amr_finder_tojoin <- amr_finder %>%
  select(c("seq_name", "amr_genes")) %>%
  mutate(seq_name = gsub(seq_name, pattern = "-", replacement = "\\."))

colnames(amr_finder_tojoin) <- c("feature", "amr_genes")



# Add prokka results
# Join with prokka output
library(ape)
prokka <- read.gff("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/prokka_result_all_clean.gff")


prokka <- prokka %>%
  mutate(Product = str_extract(attributes, pattern = "product=(.*)")) %>%
  mutate(Product = gsub(Product, pattern = "product=", replacement = ""))%>%
  mutate(Name = str_extract(attributes, pattern = "Name=(.*?);"))%>%
  mutate(Name = gsub(Name, pattern = "Name=", replacement = "")) %>%
  mutate(Name = gsub(Name, pattern = ";", replacement = "")) %>%
  mutate(Gene = str_extract(attributes, pattern = "gene=(.*?);")) %>%
  mutate(Gene  = gsub(Gene , pattern = "gene=", replacement = ""))%>%
  mutate(Gene  = gsub(Gene , pattern = ";", replacement = ""))

prokka_tojoin <- prokka %>%
  select(c("seqid", "Product", "Gene")) 

colnames(prokka_tojoin) <- c("Contig", "Product", "Gene")

prokka_resistance <- prokka_tojoin %>%
  filter(str_detect(Product, pattern = "resistance") |
           str_detect(Product, pattern = "Ars") |
           str_detect(Product, pattern = "Beta-lactamase")|
           str_detect(Product, pattern = "Silver") |
           str_detect(Product, pattern = "Zinc")  ) %>%
  mutate(Resistance = case_when(str_detect(Product, pattern = "drug") ~ "ARG",    # antibiotics
                                str_detect(Product, pattern = "Colistin") ~ "ARG",
                                str_detect(Product, pattern = "Antiseptic") ~ "ARG",
                                str_detect(Product, pattern = "Fosfomycin") ~ "ARG",
                                str_detect(Product, pattern = "Vancomycin") ~ "ARG",
                                str_detect(Product, pattern = "Daunorubicin") ~ "ARG",
                                str_detect(Product, pattern = "antibiotic") ~ "ARG",
                                str_detect(Product, pattern = "Beta") ~ "ARG",
                                str_detect(Product, pattern = "Copper") ~ "HMRG",    # metals
                                str_detect(Product, pattern = "Cobalt") ~ "HMRG",
                                str_detect(Product, pattern = "Merc") ~ "HMRG",
                                str_detect(Product, pattern = "Ars") ~ "HMRG",
                                str_detect(Product, pattern = "Nickel") ~ "HMRG",
                                str_detect(Product, pattern = "Silver") ~ "HMRG",
                                str_detect(Product, pattern = "Zinc") ~ "HMRG",
                                TRUE ~ "Other")) %>%
  mutate(feature = gsub(Contig, pattern = "-", replacement = "\\.")) 

prokka_arg <- prokka_resistance %>%
  filter(Resistance == "ARG")

prokka_hmrg <- prokka_resistance %>%
  filter(Resistance == "HMRG")



# Add kegg results
#ko2level <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/ko2level_Jan2021.csv")
ko2level <- read_csv("/Volumes/BunnyBike/mge_urban/local/kegg_ko_complete.csv")

plasmids_kegg <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plasmids_kegg_results.csv")


plasmids_kegg_level <- plasmids_kegg %>%
  mutate(KO = query_name) %>%
  left_join(., ko2level, by = "KO")

plasmids_kegg_level_best <- plasmids_kegg_level %>%
  dplyr::group_by(target_name) %>% # group by same read name
  dplyr::arrange(e_value,.by_group = T) %>% # sort by evalue within group
  slice_head(n = 1) # select the first hit (best one) for each group 

xenob_kegg <- plasmids_kegg_level_best %>%
  filter(level2 == "09111 Xenobiotics biodegradation and metabolism") %>%
  mutate(Contig = gsub(target_name, pattern = "_[0-9]$", replacement = ""))%>%
  mutate(Contig = gsub(Contig, pattern = "-", replacement = "\\."))


# Add Bacmet


ptus_bacmet <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/ptus_bacmet.out", sep = "\t", header = F)

colnames(ptus_bacmet) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                           "qend", "sstart", "send", "evalue", "bitscore")
ptus_bacmet <- ptus_bacmet %>%
  distinct(qseqid, .keep_all = T) %>%
  mutate(Contig = gsub(qseqid, pattern = "-", replacement = "\\."))



sig_families_plot <- sig_families %>%
  mutate(color_code = case_when(coef > 0 ~ "green", 
                                coef < 0 ~ "orange")) %>%
  # add abricate
  mutate(abricate = case_when(feature == "SC.T1.3_k127_822506" ~ 1, 
                          TRUE ~ 0)) %>%
  mutate(conjugation = case_when(feature %in% conjugation_tojoin$feature ~ "yes", 
                                 TRUE ~ "no"))%>%
  mutate(amr_finder = case_when(feature %in% amr_finder_tojoin$feature ~ 1, 
                                 TRUE ~ 0)) %>%
  mutate(arg_prokka = case_when(feature %in% prokka_arg$feature~ 1, 
                                TRUE ~ 0)) %>%
  mutate(hmrg_prokka = case_when(feature %in% prokka_hmrg$feature ~ 1, 
                                TRUE ~ 0)) %>%
  mutate(xenob_kegg = case_when(feature %in% xenob_kegg$Contig~ 1, 
                                TRUE ~ 0)) %>%
  mutate(bacmet = case_when(feature %in% ptus_bacmet$Contig~ 1, 
                                TRUE ~ 0)) %>%
  mutate(arg_yes = case_when(amr_finder + abricate + arg_prokka > 0 ~ 1, 
                             TRUE ~ 0)) %>%
  mutate(hmrg_yes = case_when(hmrg_prokka + bacmet > 0 ~ 1, 
                             TRUE ~ 0)) %>%
  mutate(resistance = case_when(xenob_kegg + hmrg_yes + arg_yes > 2 ~ "ARG+HMRG+XENOB", 
                                arg_yes + hmrg_yes > 1 ~ "ARG+HMRG", 
                                arg_yes + xenob_kegg > 1 ~ "ARG+XENOB", 
                                xenob_kegg + hmrg_yes > 1 ~ "HMRG+XENOB", 
                                xenob_kegg > 0 ~ "XENOB", 
                                hmrg_yes > 0 ~ "HMRG", 
                                arg_yes > 0 ~ "ARG",
                                TRUE ~ "None"))  %>%
  mutate(number_resistance_genes = amr_finder + abricate + arg_prokka + hmrg_prokka + xenob_kegg + bacmet) 

library(see)
library(ggpattern)

sig_families_plot$resistance <- factor(sig_families_plot$resistance, levels = c("ARG", "HMRG","XENOB", "ARG+HMRG","HMRG+XENOB","ARG+HMRG+XENOB", "None"))
sig_families_plot$resistance <- factor(sig_families_plot$resistance, levels = c("ARG", "HMRG","XENOB", "ARG+HMRG","HMRG+XENOB", "None"))
sig_families_plot$conjugation <- as.factor(sig_families_plot$conjugation)

ggplot(sig_families_plot, aes(x = reorder(feature, -coef), y = coef, fill = resistance, pattern = conjugation))+
  #geom_col() +
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.010,
                   pattern_key_scale_factor = 0.06) +
  geom_text(aes(label = number_resistance_genes), vjust = -0.5, size = 2) +
  scale_fill_manual(name = "Resistance", values = c("deeppink1","orange", "forestgreen" ,"purple", "deepskyblue","chocolate4", "grey86"))+
  #scale_fill_see()+
  scale_pattern_manual(name = "ConjG", values = c(yes = "stripe", no = "none")) +
  xlab("Plasmid") +
  ylab("Coefficient") +
  theme_bw() +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+
  theme(legend.position = "right", axis.text.x=element_blank())
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmids_maaslin_resistance_conjugation_mixed.pdf", device = "pdf", width = 9 , height = 5, units = "in")



# ------------------------------------------------------------------------------------------------------------------------------
## Mantel tests

# To run this you have to load the distance matrices for each thing on their script

mantel_result <- mantel(ptus.bray, bracken.bray, method = "pearson", permutations = 999)


#  ------------------------------------------------------------------------------------------------------------------------------
# testing the effect of landuse

pal2 <- c("#ef3e36","#17bebb","#2e282a","#fad8d6")
levels(as.factor(md_final$landuse2))

# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("landuse2", "ptu_shannon")) %>%
  na.omit() %>%
  dplyr::group_by(landuse2) %>%
  dplyr::summarize(sd = sd(ptu_shannon),
                   ptu_shannon  = mean(ptu_shannon))

ggplot(md_final, aes(x = landuse2, y = ptu_shannon)) + 
  geom_jitter(aes(fill= landuse2, color = landuse2), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin = ptu_shannon-sd, ymax = ptu_shannon+sd, color = landuse2),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal2)+
  scale_color_manual(values=pal2)+
  xlab(NULL) + 
  ylab("PTU Shannon diversity")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptu_shannon_ecosystem.pdf", device = "pdf", width = 5, height = 5 , units = "in")

md_nourban <- md_final %>%
  filter(!urban.natural == "urban")

model <- lmer(ptu_rich ~ landuse2 + (1|site), data = md_nourban)
model <- lmer(ptu_shannon ~ landuse2 + (1|site), data = md_nourban)
model <- lmer(ptus_relab ~ landuse2 + (1|site), data = md_nourban)
Anova(model, type = 3)




# PERMANOVA

# select samples only in non urban
md_nourban <- md_final %>%
  filter(!urban.natural == "urban")

# select samples on ptu table non urban

rpkm_ptus_final_nourban <- rpkm_ptus_final %>%
  select(md_nourban$sample)

ptus.bray_nourban <- vegdist(t(rpkm_ptus_final_nourban), method="bray")

#ordination (non-multidimensional scaling)
ptus.nmds <- metaMDS(ptus.bray_nourban, k=2, try = 100)
md_nourban$Axis01 = ptus.nmds$points[,1]
md_nourban$Axis02 = ptus.nmds$points[,2]
ptus.nmds$stress #0.1014412

pal3 <- c("#ef3e36","#17bebb","#2e282a")
# urban.natural
ggplot(md_nourban, aes(Axis01, Axis02))+
  geom_point(aes(color=landuse2), size=4)+
  #stat_ellipse(aes(color = landuse2))+
  scale_color_manual(values=pal3)+
  theme_bw()+
  theme(legend.position="bottom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_beta_ecosystem_nourban.pdf", device = "pdf", width = 5, height = 5 , units = "in")


adonis2(ptus.bray_nourban   ~ landuse2, data = md_nourban, permutations = 999, method = "bray")


## dispersion
# Calculate beta dispersion
ptus.disp <- betadisper(ptus.bray_nourban, md_nourban$landuse2)
# Test for significant differences in dispersion
permutest(ptus.disp, permutations = 999)
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
  scale_fill_manual(values = pal3) +
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
  scale_color_manual(values = pal3) +
  theme_bw() +
  labs(title = "PCoA of Dispersion")

p2 + p1 + plot_layout(widths = c(1.2, 1))

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/ptus_comp_dispersion_nourban.pdf", 
       width = 10, height = 4, units = "in")


# Anova for testing
anova(ptus.disp)


