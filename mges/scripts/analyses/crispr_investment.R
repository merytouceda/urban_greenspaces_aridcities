# CRISPR investment analyses for urban plasmids project

library(lme4)
library(tidyverse)
library(car)
library(ggpubr)
pal1 = c("#FF823B", "#359951")


########################
#### Calculate CRIPR Investment
########################
# CRISPR investment =  number of CRISPR arrays / median number of copies of ribosomal proteins found in a sample.

# average copy number
acn <- read.table("/Volumes/BunnyBike/landuse_functionaltraits/data/sample_16S_coverage.txt", header = T)

acn <- acn %>%
  rownames_to_column(var = "sample") %>%
  mutate(sample = gsub(sample, pattern = "-", replacement = "\\."))

colnames(acn) <- c("sample", "acn")

# CRISPR number in md
md <- read_csv("/Volumes/BunnyBike/landuse_functionaltraits/data/LU_metadata.csv")

# Calculate investment
md_crispr <- md %>%
  left_join(., acn, by = "sample") %>%
  mutate(crispr_i = CRISPR/acn)


# plot
# Line with error bar
mrelab_lineplot <- md_crispr %>%
  dplyr::select(c("urban.natural", "crispr_i")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural) %>%
  dplyr::summarize(sd = sd(crispr_i),
                   crispr_i = mean(crispr_i))

ggplot(md_crispr, aes(x = urban.natural, y = crispr_i)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin =  crispr_i-sd, ymax =  crispr_i+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("CRISPR Investment")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/crispr_investment.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(crispr_i ~ urban.natural + (1|site), data = md_crispr), type = 3)
#Response: crispr_i
#Chisq Df Pr(>Chisq)    
#(Intercept)   22.9952  1  1.624e-06 ***
#  urban.natural  0.3307  1     0.5652  


########################
#### Correlation with AMR abundance
########################
# Run the amrfinder_abricate.R script to get the abundance of AMRFinder and Abricate AMR genes

# join the mds
md_crispr_join <- md_crispr %>%
  mutate(sample = gsub(sample, pattern = "\\.", replacement = "-")) %>%
  select(c("sample", "crispr_i"))

md_final_crispr <- md_final %>%
  left_join(., md_crispr_join, by = "sample")

# correlation
cor.test(md_final_crispr$ptus_relab, md_final_crispr$crispr_i)

ggplot(md_final_crispr, aes(x = ptus_relab, y = crispr_i))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", linetype = "dashed", color = "grey30", fill = "grey80")+
  stat_cor(method = "pearson")+
  xlab("Relative abundance PTUs") + 
  ylab("CRISPR Investment")+
  theme_bw() #+
  facet_wrap(~urban.natural)
  ggsave("/Volumes/BunnyBike/mge_urban/local/figures/crispr1_ptusrelab_corr.pdf", device = "pdf", width = 5, height = 5 , units = "in")
  

  
  ########################
  #### Correlation with plasmid abundance
  ########################
  # join the mds
  md_crispr_join <- md_crispr %>%
    select(c("sample", "crispr_i"))
  
  md_final_crispr <- md_final %>%
    left_join(., md_crispr_join, by = "sample")
  
  # correlation
  cor.test(md_final_crispr$ptus_relab, md_final_crispr$crispr_i)
  
  ggplot(md_final_crispr, aes(x = amr_relab, y = crispr_i))+
    geom_point(size = 2)+
    geom_smooth(method = "lm", linetype = "dashed", color = "grey30", fill = "grey80")+
    stat_cor(method = "pearson")+
    xlab("Relative abundance of AMR genes") + 
    ylab("CRISPR Investment")+
    theme_bw() #+
  facet_wrap(~urban.natural)
  ggsave("/Volumes/BunnyBike/mge_urban/local/figures/crispr1_amrrelab_corr.pdf", device = "pdf", width = 5, height = 5 , units = "in")
  
  ########################
  #### Correlation with MGE abundance
  ########################
  # join the mds
  md_crispr_join <- md_crispr %>%
    select(c("sample", "crispr_i"))
  
  md_final_crispr <- md %>%
    left_join(., md_crispr_join, by = "sample")
  
  # correlation
  cor.test(md_final_crispr$mge_relab, md_final_crispr$crispr_i)
  
  ggplot(md_final_crispr, aes(x = mge_relab, y = crispr_i))+
    geom_point(size = 2)+
    geom_smooth(method = "lm", linetype = "dashed", color = "grey30", fill = "grey80")+
    stat_cor(method = "pearson")+
    xlab("Relative abundance of MGTUs") + 
    ylab("CRISPR Investment")+
    theme_bw() #+
  facet_wrap(~urban.natural)
  ggsave("/Volumes/BunnyBike/mge_urban/local/figures/crispr1_mgturelab_corr.pdf", device = "pdf", width = 5, height = 5 , units = "in")
  