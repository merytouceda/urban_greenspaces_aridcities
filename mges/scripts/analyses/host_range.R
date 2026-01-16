# analysis of the Host range data of urban plasmids

library(tidyverse)
library(vegan)
library(car)
library(performance)
library(lme4)

pal1 = c("#FF823B", "#359951")
# Load and clean 
host_range_s <- read_tsv("/Volumes/BunnyBike/mge_urban/local/genomad_on_all_contigs/HostRange_Species.tsv")

separations <- as.character(c(1:205))

host_range_s_clean <- host_range_s %>%
  add_row(!!!setNames(colnames(.), colnames(.)), .before = 1) %>%
  # Give new generic column names
  set_names(c("plasmid", paste0("col", 2:ncol(.)))) %>%
  # Remove any asterisks from the first row values
  mutate(across(everything(), ~if_else(row_number() == 1, 
                                       str_remove_all(., "\\*"), 
                                       .))) %>%
  mutate(col2 = gsub(col2, pattern = "f__", replacement = "\\| f__")) %>%
  mutate(col2= gsub(col2, pattern = "^\\|", replacement = "")) %>%
  filter(!col2 == "No species was predicted") %>%
  mutate(degree = case_when(str_count(col2, pattern = "f__") < 2 ~ 'I, II, III', 
                           TRUE ~ 'IV, V, VI'))


# % of plasmids with broad host range
percent_range <- host_range_s_clean %>%
  group_by(degree) %>%
  summarize(n = n())



##### PLOT PRESENCE OF DEGREE IN URBAN AND NATURAL

# Join this information with presence absence in urban and natural

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

host_range_count <- host_range_s_clean %>%
  select(c("plasmid","degree")) %>%
  left_join(., ptu_pab, by = "plasmid") %>%
  drop_na() %>%
  select(-plasmid) %>%
  pivot_longer(-degree, names_to = "sample", values_to="count") %>%
  group_by(degree, sample) %>%
  summarize(counts = sum(count)) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural")) %>%
  mutate(site = case_when(str_detect(sample, pattern = "RP") ~ "RP", 
                          str_detect(sample, pattern = "HP") ~ "HP",
                          str_detect(sample, pattern = "SC") ~ "SC",
                          str_detect(sample, pattern = "RC") ~ "RC",
                          str_detect(sample, pattern = "RP") ~ "RP",
                          str_detect(sample, pattern = "11") ~ "11",
                          str_detect(sample, pattern = "8") ~ "8",
                          str_detect(sample, pattern = "UAB") ~ "UAB",
                          str_detect(sample, pattern = "Ex") ~ "Ex",

))

# for error line 
# Line with error bar
mrelab_lineplot <- host_range_count  %>%
  ungroup() %>%
  dplyr::select(c("urban.natural", "counts", "degree")) %>%
  na.omit() %>%
  dplyr::group_by(urban.natural, degree) %>%
  dplyr::summarize(sd = sd(counts),
                   counts = mean(counts))

ggplot(host_range_count, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~degree)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number plasmids")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/host_range_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")


ggplot(host_range_count, aes(x = urban.natural, y = counts)) + 
  geom_jitter(aes(fill= urban.natural, color = urban.natural), size = 2.5, alpha = 0.6)+
  #geom_errorbar(aes(ymin = mrelab_lineplot$av_relab-mrelab_lineplot$sd, ymax = mrelab_lineplot$av_relab+mrelab_lineplot$sd), width=.1,
  #position=position_dodge(0.05), color = "black", linewidth = 0.4, linetype = "dashed")+
  geom_pointrange(aes(ymin = counts-sd, ymax = counts+sd, color = urban.natural),size = 0.5, linewidth = 1, data = mrelab_lineplot)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number of plasmids")+
  theme_bw()+
  facet_wrap(~degree)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/host_range_pointrange.pdf", device = "pdf", width = 7, height = 5 , units = "in")


hr_narrow <- host_range_count %>%
  filter(degree == "I, II, III")
Anova(lmer(counts ~ urban.natural + (1|site), data = hr_narrow ), type = 3)

hr_broad <- host_range_count %>%
  filter(!degree == "I, II, III")

Anova(lmer(counts ~ urban.natural + (1|site), data = hr_broad), type = 3)

Anova(lm(counts ~ degree * urban.natural, data = host_range_count))



# Now plot the abundance of those plasmids
# NOTE: RUN THE RPKM CALCULATION FROM ptus_analyses.R
ptu_counts_join <- rpkm_ptus_final %>%
  rownames_to_column(var = "plasmid") %>%
  pivot_longer(-c("plasmid"), names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "-", replacement = "\\.")) %>%
  pivot_wider(names_from = "sample", values_from = "count")

host_range_rpkm <- host_range_s_clean %>%
  select(c("plasmid","degree")) %>%
  left_join(., ptu_counts_join, by = "plasmid") %>%
  drop_na() %>%
  select(-plasmid) %>%
  pivot_longer(-degree, names_to = "sample", values_to="count") %>%
  group_by(degree, sample) %>%
  summarize(counts = sum(count)) %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural"))

ggplot(host_range_rpkm, aes(x = urban.natural, y = counts))+
  geom_boxplot(aes(color = urban.natural))+
  geom_jitter(aes(color = urban.natural))+
  facet_wrap(~degree)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Abundance (RPKM)")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/host_range_rpkm_urbannatural.pdf", device = "pdf", width = 8, height = 8 , units = "in")







