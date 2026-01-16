# Unknown kegg analysis of plasmids from urban project


library(tidyverse)
library(see)
library(pals)
library(vegan)
library(Polychrome)
library(car)
library(lme4)
library(ggExtra)

pal1 = c("#56282D", "#77966D")

# load metadata
md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")

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


# ORF counts
orf_counts <- read.table("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plasmid_genes_count_table.txt", sep = "\t", header = T)

unk_count <- orf_counts %>%
  filter(!Contig %in% plasmids_kegg_level_best$target_name) %>%
  column_to_rownames(var = "Contig")

# Number of unknown genes = 1397 

#####
## Richness
#### 
summary(colSums(unk_count))

md$unk_func_rich <- specnumber(rrarefy(t(unk_count), sample = 634))

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "unk_func_rich")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(unk_func_rich),
                   unk_func_rich  = mean(unk_func_rich))

p <- ggplot(md, aes(x = urban.natural, y = unk_func_rich)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  unk_func_rich-sd, ymax =  unk_func_rich+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Plasmid unknown functional richness")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_unk_functional_richness_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(unk_func_rich ~ urban.natural + (1|site), data = md), type = 3)
#Chisq Df Pr(>Chisq)    
#(Intercept)   216.905  1     <2e-16 ***
#  urban.natural   0.591  1      0.442  


#######################
# RPKM
#######################
# length
len <- read_csv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/plasmid_genes_urban_rep_len.csv")
colnames(len) <- c("Contig", "len")

# select the total reads per sample from the metadata
total_reads <- md %>%
  dplyr::select(sample, total.reads)

unk_count_long <- unk_count %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "Count") %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = ""))

# 3. join the total reads to each row based on sample
unk_count_long_plusreads <- left_join(unk_count_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
unk_count_long_plusreads_plustlen <- left_join(unk_count_long_plusreads, len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_unk_final <- unk_count_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(total.reads*as.numeric(len))) %>%
  dplyr::select(-c("Count", "total.reads", "len")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) %>%
  column_to_rownames(var = "Contig")



############################
# beta-diversity
############################

kegg.bray <- vegdist(t(rpkm_unk_final), method="bray")

#ordination (non-multidimensional scaling)
mges.nmds <- metaMDS(kegg.bray, k=2, try = 100)
md$Axis01 = mges.nmds$points[,1]
md$Axis02 = mges.nmds$points[,2]
mges.nmds$stress #0.03989818 the smaller the better, good <0.3)

ggplot(md, aes(Axis01, Axis02))+
  geom_point(aes(color=urban.natural), size=4)+
  stat_ellipse(aes(color = urban.natural))+
  scale_color_manual(values=pal1)+
  theme_bw()+
  theme(legend.position="botoom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmids_unk_profile.pdf", device = "pdf", width = 5, height = 5 , units = "in")

adonis2(kegg.bray  ~ urban.natural, data = md, permutations = 999, method = "bray", strata = md$site)
adonis2(kegg.bray  ~ urban.natural, data =  md, permutations = 999, method = "bray") 
#Df SumOfSqs      R2      F Pr(>F)    
#urban.natural  1    1.514 0.18868 5.1163  0.001 ***



############################
# abundance
############################
sum_rel_ab_unk <- as.data.frame(t(rpkm_unk_final)) 

sum_rel_ab_unk_final <- sum_rel_ab_unk  %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(unk_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(unk_relab) %>%
  rownames_to_column(var = "sample")

md_final <- md %>%
  left_join(., sum_rel_ab_unk_final, by = "sample")


# Line with error bar
mrelab_lineplot <- md_final %>%
  dplyr::select(c("urban.natural", "unk_relab")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(unk_relab),
                   unk_relab  = mean(unk_relab))

ggplot(md_final, aes(x = urban.natural, y = unk_relab)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin = unk_relab-sd, ymax = unk_relab+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Unknown gene abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/plasmid_unk_relab_pointrange.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(unk_relab ~ urban.natural + (1|site), data = md_final), type = 3)
#Chisq Df Pr(>Chisq)   
#(Intercept)   4.2718  1   0.038750 * 
#  urban.natural 9.3653  1   0.002211 **

Anova(lm(unk_relab ~ urban.natural * vegetation_structure, data = md_final))


################################################################################################################
# using ptu counts ? Ask alise (get from kegg_plasmid_analysis)
#########################
#Alise: Calculate ratio of unknowns per plasmid, plot the differences in ratio urban vs. natural


# make a list of the unknowns

unk_genes <- unk_count %>%
  rownames_to_column(var = "gene") %>%
  select(gene) %>%
  mutate(plasmid = gsub(gene, pattern = "_[0-9]*$", replacement = ""))

unk_per_plasmid <- unk_genes %>%
  group_by(plasmid) %>%
  summarize(unk_counts = n())


# now create a table with the number of known genes per plasmid
known_per_plasmid <- plasmids_kegg_level_best %>%
  select(target_name) %>%
  mutate(plasmid = gsub(target_name, pattern = "_[0-9]*$", replacement = "")) %>%
  group_by(plasmid) %>%
  summarize(known_counts = n())

# merge 
plasmid_unk_known <- unk_per_plasmid %>%
  left_join(., known_per_plasmid, by = "plasmid") %>%
  replace(is.na(.), 0) %>%
  mutate(ratio_unk = unk_counts/(unk_counts + known_counts))

summary(plasmid_unk_known$ratio_unk)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.04545 0.27273 0.44444 0.43734 0.58333 1.00000 
hist(plasmid_unk_known$ratio_unk)

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

# mean unk ratio per sample
unk_ratio_sample <- ptu_pab  %>%
  pivot_longer(-plasmid, names_to = "sample", values_to = "pab") %>%
  left_join(., plasmid_unk_known, by = "plasmid") %>%
  filter(!pab == 0) %>%
  drop_na() %>%
  group_by(sample) %>%
  summarize(mean_unkratio = mean(ratio_unk))

mean(md$mean_unkratio)

md <- md %>%
  left_join(unk_ratio_sample, by = "sample")

# Line with error bar
mrelab_lineplot <- md %>%
  dplyr::select(c("urban.natural", "mean_unkratio")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(mean_unkratio),
                   mean_unkratio = mean(mean_unkratio))

p <- ggplot(md, aes(x = urban.natural, y = mean_unkratio)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  mean_unkratio-sd, ymax =  mean_unkratio+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average unknown ratio")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))

# add marginal
p_marg <- ggMarginal(p, 
                     type = "histogram",
                     margins = "y",
                     size = 8,
                     fill = "gray", color = "gray60")

p_marg

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/unk_ratio_marg.pdf", p_marg, device = "pdf", width = 5, height = 5 , units = "in")


Anova(lmer(mean_unkratio ~ urban.natural + (1|site), data = md), type = 3)
#Chisq Df Pr(>Chisq)    
#(Intercept)   11728.3085  1     <2e-16 ***
#  urban.natural     0.1993  1     0.6553






# --------------------
## Model diagnostics
# models
model <- lmer(unk_func_rich ~ urban.natural + (1|site), data = md)
model <- lmer(unk_relab ~ urban.natural + (1|site), data = md_final)
model <- lmer(mean_unkratio ~ urban.natural + (1|site), data = md)
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

ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mean_unkratio_model_diagnostics.pdf", 
       width = 10, height = 4, units = "in")

# Shapiro-Wilk test for normality
shapiro.test(residuals(model))
# unk rich: W = 0.96369, p-value = 0.5169
# unk relab: W = 0.77896, p-value = 0.000136
# mean unk ratio: W = 0.98158, p-value = 0.9231

# Check model (performance)
check_model(model, check=c("linearity", "homogeneity", "normality"))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mean_unkratio_checkmodel.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# DO PERMUTATION TEST unk relab

# permutation test: 
model_perm <- lmp(unk_relab ~ urban.natural, data = md_final)
summary(model_perm)
# urban.natural1   -3.039 5000   0.0056 **
