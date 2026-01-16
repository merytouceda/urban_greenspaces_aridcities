# BARPLOT TYPES OF MGES

library(tidyverse)
library(vegan)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(see)
pal1 = c("#FF823B", "#359951")

library(RColorBrewer)

pal2 <- brewer.pal(4, "Set3")

##############
## BLAST OUTPUT FILTER
##############

mobileog_diamond <- read.table("/Volumes/BunnyBike/mge_urban/local/mobileOG/local_urban_mobileOG_fixednames_100124.out", header = F)

colnames(mobileog_diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                                "qend", "sstart", "send", "evalue", "bitscore")


# keep only the best hit per gene
mobileog_diamond_filtered <- mobileog_diamond %>%
  dplyr::group_by(qseqid) %>% # group by same read name
  dplyr::arrange(bitscore,.by_group = T) %>% # sort by evalue within group
  slice_head(n = 1) %>% # select the first hit (best one) for each group 
  filter(!pident < 90) %>%
  filter(length > 50)

mean(mobileog_diamond_filtered$length)
summary(mobileog_diamond_filtered$length)
hist(mobileog_diamond_filtered$length)

final_diamond <- mobileog_diamond_filtered  %>%
  separate(sseqid, into = c("mobileOG_id", "gene_name", "other_id", "category", "something", "db", "class"), sep = "\\|") %>%
  mutate(MGE_type = case_when(db %in% c("ISFinder") ~ "Insertion_sequence", 
                              db %in% c("AICE","ICE","CIME","IME","immedb") ~ "Integrative_element",
                              db %in% c("COMPASS","Plasmid") ~ "Plasmid",
                              db %in% c("ACLAME", "Multiple") ~ "Multiple",
                              db %in% c("pVOG","GPD") ~ "Phage")) %>%
  filter(!MGE_type == "Phage") %>%
  filter(!category == "phage") 

final_diamond_plot <- final_diamond %>%
  group_by(MGE_type, category) %>%
  summarize(count = n()) %>%
  replace_na(list(MGE_type = "Other"))


ggplot(final_diamond_plot, aes(x = reorder(category, count), y = count, fill = MGE_type))+
  geom_col()+
  xlab("MGE-associated molecular machinery") + 
  ylab("Number of genes")+
  coord_flip()+
  scale_fill_manual(values = pal2)+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mge_types.pdf", device = "pdf", width = 10, height = 4 , units = "in")


final_diamond_pie <- final_diamond_plot %>%
  group_by(MGE_type) %>%
  summarise(
    total_count = sum(count),
    percentage = round(total_count / sum(final_diamond_plot$count) * 100, 2)
  ) %>%
  mutate(
  # Add percentage sign
  label = paste0(round(percentage, 1), "%"),
  # Calculate the position for labels
  position = cumsum(percentage) - percentage/2, 
  pos_rad = position * pi / 180
)
            
ggplot(final_diamond_pie, aes(x = 2, y = percentage, fill = MGE_type)) +
  # Create the donut
  geom_bar(stat = "identity", width = 1) +
  # Add percentage labels
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), color = "White") +
  geom_text(aes(label = MGE_type),
            position = position_stack(vjust = 0.4), color = "White") +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y", start = 0) +
  # Remove unnecessary elements and customize theme
  theme_void() +
  scale_fill_manual(values = pal2)+
  # Create a hole in the center to make it a donut
  xlim(0.5, 2.7) +
  # Custom colors and legend position
  theme(legend.position = "bottom", plot.margin = margin(1, 1, 1, 1, "cm")) 
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mge_types_pie.pdf", device = "pdf", width = 5, height = 5 , units = "in")

  


## Do the same plots but with abundance

# Run the RPKM calculation from mgtus_analyses

# add contig to diamond result
final_diamond_tojoin <- final_diamond %>%
  mutate(MGE_type = str_trim(MGE_type)) %>%
  mutate(Contig = gsub(qseqid, pattern = "_[0-9]*$", replacement = "")) %>%
  mutate(Contig = gsub(Contig, pattern = "\\./", replacement = "")) %>%
  ungroup() %>%
  select(c("Contig", "category", "MGE_type")) %>%
  group_by(Contig) %>%
  summarize(number_mges = n(),
              category = if(n_distinct(category) == 1) first(category) else "Multiple",
            MGE_type = {
              types <- unique(MGE_type[MGE_type != "Multiple"])
              if(length(types) == 1) types else "Multiple"
            })

check <- final_diamond %>%
  mutate(Contig = gsub(qseqid, pattern = "_[0-9]*$", replacement = "")) %>%
  mutate(Contig = gsub(Contig, pattern = "\\./", replacement = "")) %>%
  group_by(Contig) %>%
  summarise(n_distinct_types = n_distinct(MGE_type),
            all_types = paste(unique(MGE_type), collapse = ", ")) %>%
  filter(n_distinct_types > 1)



# join to counts
mgtus_type_counts <- rpkm_mgtus_final  %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "count") %>%
  left_join(.,final_diamond_tojoin, by = "Contig") %>%
  mutate(urban.natural = case_when(sample %in% c("HP.T1.3",   "HP.T1.6", "HP.T1.9" ,
                                                 "RP.T1.3" ,  "RP.T1.6",   "RP.T1.9") ~ "urban", 
                                   TRUE ~ "natural")) %>%
  group_by(category, urban.natural, MGE_type) %>%
  summarize(counts = sum(count))
 

ggplot(mgtus_type_counts, aes(x = reorder(category, -counts), y = counts, fill = MGE_type))+
  geom_col()+
  xlab("MGE-associated molecular machinery") + 
  ylab("Abundance (RPKM)")+
  scale_fill_manual(values= pal2)+
  theme_classic()+
  #facet_wrap(~urban.natural)+
  theme(legend.position = "right", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mge_types_abundance_urbannatural.pdf", device = "pdf", width = 10, height = 4 , units = "in")


# NOTE: maybe repeat these with the abundance of the gene? 
