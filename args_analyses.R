# ARGs analyses, land use project
# Mery Touceda-Suárez
# September 2022


# load libraries
library(tidyverse)
library(ggplot2)
library(stringr)
library(dplyr)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)
#install.packages("VennDiagram")
library(VennDiagram)

setwd("/Volumes/MyPassport/githubrepos/LandUse")

# load and process data
part1 <- read.table("card_part1.out.txt", sep ="\t", header = T, fill = T)
part2 <- read.table("card_part2.out.txt", sep ="\t", header = T, fill = T)
part3 <- read.table("card_part3.out.txt", sep ="\t", header = T, fill = T)
metadata <- read.table("metadata-landuse.txt", sep ="\t", header = T)

# concatenate all parts
allparts <- rbind(part1, part2, part3)

#remove empty columns
allparts<- allparts[,-c(2,3,4,5)]

#take sample name
allparts$sample<-str_extract(allparts$ORF_ID, "[^_]+")

#take ORF number
## problem! variable lengh on ORF number    # SOLVE THIS
#part1$ORF_number = substr(part1$ORF_ID,9,)

#get rid of some spurious rows
allparts<- allparts[!allparts$sample=="100.0",]
allparts<- allparts[!allparts$sample=="n/a",]
allparts<- allparts[!allparts$sample=="3004680",]
allparts<- allparts[!allparts$sample=="protein homolog model",]
allparts<- allparts[!allparts$sample=="aminoglycoside antibiotic",]
allparts<- allparts[!allparts$sample=="3002593",]
allparts<- allparts[!allparts$sample=="antibiotic inactivation",]
allparts<- allparts[!allparts$sample=="AAC(2)\t\tNALRCGYLLVGAAERVYVRVTGDEVARV
                    PRYEDDAVHQLLRRRWLTTGASHSLTCGAASLTGIALLVPKQTRAAVARWDFLQRPPSWPQQRGTT
                    TTEPGRPDQDPRAGRVVRLDDHRRRRR*\tMDTHHVHTARLVHTADLDGETLRRLQQMVTDAFAGDFD
                    ETDWEHALGGMHALIWRHGTIIAHAAVVQRRLFYHGNALRCGYLEGVAVRKDCRGRGLVHALLDAIEQVI
                    RGAYQFGALSSSDRARRVYMSRGWLPWLGPTSVLAPTGVIRTPDDDGSVFVLPVGINPDTSSGLMCDWRAG
                    NVW\t67.03\tgnl|BL",]



#count number of different ARGs per sample
number_args = c()
number_args = c(number_args, length(unique(allparts[allparts$sample=="11B-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="11B-T1-6",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="11B-T1-9",]$Best_Hit_ARO)))

number_args = c(number_args, length(unique(allparts[allparts$sample=="8-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="8-T1-6",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="8-T1-9",]$Best_Hit_ARO)))

number_args = c(number_args, length(unique(allparts[allparts$sample=="Ex45-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="Ex45-T1-6",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="Ex45-T1-9",]$Best_Hit_ARO)))

number_args = c(number_args, length(unique(allparts[allparts$sample=="HP-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="HP-T1-6",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="HP-T1-9",]$Best_Hit_ARO)))

number_args = c(number_args, length(unique(allparts[allparts$sample=="RC-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, NA) #not present in card results table
number_args = c(number_args, NA) #not present in card results table

number_args = c(number_args, NA) #not present in card results table
number_args = c(number_args, length(unique(allparts[allparts$sample=="RP-T1-6",]$Best_Hit_ARO)))
number_args = c(number_args, length(unique(allparts[allparts$sample=="RP-T1-9",]$Best_Hit_ARO)))

number_args = c(number_args, length(unique(allparts[allparts$sample=="SC-T1-3",]$Best_Hit_ARO)))
number_args = c(number_args, NA) #not present in card results table
number_args = c(number_args, NA) #not present in card results table

number_args = c(number_args, NA) #not present in card results table
number_args = c(number_args, NA) #not present in card results table
number_args = c(number_args, NA) #not present in card results table

#add as a column to metadata to plot
metadata$args_number <- number_args

#quick box plot
ggplot(metadata, aes(x = land.use, y = args_number, fill = land.use))+
  geom_boxplot(alpha=0.8)+
  scale_fill_manual(values = c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) +
  ylab("Number of different ARGs")+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("arg_rich.pdf", device = "pdf", width = 7, height = 5 , units = "in")



# next steps: 
# QUESTION : should i do the unique value of each sample, or should I put them all together by land use and then do the unique value?
# 1. match the ARG to the "abundance" of such contig or ORF, to see the most abundant ARGs in each sample


# 2. take the list of ARGs in the different land use types and search their potential origin

# build dataframe to export


# PASTURE RECENT
pasture_recent <- allparts[allparts$sample %in% c("11B-T1-3","11B-T1-6","11B-T1-9",
                                                  "8-T1-3","8-T1-6","8-T1-9"),]
length(unique(pasture_recent$Best_Hit_ARO)) # number of different ARGs
unique(pasture_recent$Best_Hit_ARO)
#write.table(unique(pasture_recent$Best_Hit_ARO), "pasture_recent_args_list.txt", sep = "\t", quote=F)
unique(pasture_recent$Drug.Class)
unique(pasture_recent$Resistance.Mechanism)
unique(pasture_recent$AMR.Gene.Family)

#summarize this table by gene family
pasture_recent_AMRs <- pasture_recent[,c("Best_Hit_ARO", "Drug.Class","Resistance.Mechanism", "AMR.Gene.Family")]
pasture_recent_AMRs_grouped <- pasture_recent_AMRs %>% 
  group_by(AMR.Gene.Family) %>% 
  mutate(AMRs.in.family = paste0(Best_Hit_ARO, collapse = ",")) %>% 
  mutate(Resistance.Mechanism.in.family = paste0(Resistance.Mechanism, collapse = ",")) %>% 
  mutate(Drug.Class.in.family = paste0(Drug.Class, collapse = ",")) 
# collapse to have one row per gene family
pasture_recent_AMRs_collapsed <- pasture_recent_AMRs_grouped[,-c(1,2,3)]
pasture_recent_AMRs_collapsed <- unique(pasture_recent_AMRs_collapsed)
write.table(pasture_recent_AMRs_collapsed , "pasture_recent_arg_families.txt", sep = "\t", quote=F)


# PASTURE OLD
pasture_old <- allparts[allparts$sample %in% c("Ex45-T1-3","Ex45-T1-6","Ex45-T1-9",
                                                  "UAB-T1-3","UAB-T1-6","UAB-T1-9"),]
length(unique(pasture_old$Best_Hit_ARO))
unique(pasture_old$Drug.Class)
unique(pasture_old$Resistance.Mechanism)
unique(pasture_old$AMR.Gene.Family)
#summarize this table by gene family
pasture_old_AMRs <- pasture_old[,c("Best_Hit_ARO", "Drug.Class","Resistance.Mechanism", "AMR.Gene.Family")]
pasture_old_AMRs_grouped <- pasture_old_AMRs %>% 
  group_by(AMR.Gene.Family) %>% 
  mutate(AMRs.in.family = paste0(Best_Hit_ARO, collapse = ",")) %>% 
  mutate(Resistance.Mechanism.in.family = paste0(Resistance.Mechanism, collapse = ",")) %>% 
  mutate(Drug.Class.in.family = paste0(Drug.Class, collapse = ",")) 
# collapse to have one row per gene family
pasture_old_AMRs_collapsed <- pasture_old_AMRs_grouped[,-c(1,2,3)]
pasture_old_AMRs_collapsed <- unique(pasture_old_AMRs_collapsed)
write.table(pasture_old_AMRs_collapsed , "pasture_old_arg_families.txt", sep = "\t", quote=F)


# URBAN
urban <- allparts[allparts$sample %in% c("HP-T1-3","HP-T1-6","HP-T1-9",
                                                  "RP-T1-3","RP-T1-6","RP-T1-9"),]
length(unique(urban$Best_Hit_ARO))
unique(urban$Drug.Class)#different drug classes
length(unique(urban$Drug.Class))#number of different drug classes
unique(urban$Resistance.Mechanism)
length(unique(urban$Resistance.Mechanism))
unique(urban$AMR.Gene.Family)#different families 
length(unique(urban$AMR.Gene.Family))#number of different families
#summarize this table by gene family
urban_AMRs <- urban[,c("Best_Hit_ARO", "Drug.Class","Resistance.Mechanism", "AMR.Gene.Family")]
urban_AMRs_grouped <- urban_AMRs %>% 
  group_by(AMR.Gene.Family) %>% 
  mutate(AMRs.in.family = paste0(Best_Hit_ARO, collapse = ",")) %>% 
  mutate(Resistance.Mechanism.in.family = paste0(Resistance.Mechanism, collapse = ",")) %>% 
  mutate(Drug.Class.in.family = paste0(Drug.Class, collapse = ",")) 
# collapse to have one row per gene family
urban_AMRs_collapsed <- urban_AMRs_grouped[,-c(1,2,3)]
urban_AMRs_collapsed <- unique(urban_AMRs_collapsed)
write.table(urban_AMRs_collapsed , "urban_arg_families.txt", sep = "\t", quote=F)



# NATURAL
natural <- allparts[allparts$sample %in% c("RC-T1-3","RC-T1-6","RC-T1-9",
                                                  "SC-T1-3","SC-T1-6","SC-T1-9"),]
length(unique(natural$Best_Hit_ARO))
unique(natural$Drug.Class)#different drug classes
length(unique(natural$Drug.Class))#number of different drug classes
unique(natural$Resistance.Mechanism)
unique(natural$AMR.Gene.Family)#different families 
length(unique(natural$AMR.Gene.Family))#number of different families
#summarize this table by gene family
natural_AMRs <- natural[,c("Best_Hit_ARO", "Drug.Class","Resistance.Mechanism", "AMR.Gene.Family")]
natural_AMRs_grouped <- natural_AMRs %>% 
  group_by(AMR.Gene.Family) %>% 
  mutate(AMRs.in.family = paste0(Best_Hit_ARO, collapse = ",")) %>% 
  mutate(Resistance.Mechanism.in.family = paste0(Resistance.Mechanism, collapse = ",")) %>% 
  mutate(Drug.Class.in.family = paste0(Drug.Class, collapse = ",")) 
# collapse to have one row per gene family
natural_AMRs_collapsed <- natural_AMRs_grouped[,-c(1,2,3)]
natural_AMRs_collapsed <- unique(natural_AMRs_collapsed)
write.table(natural_AMRs_collapsed , "natural_arg_families.txt", sep = "\t", quote=F)



#Venn Diagram
# gene families
v <- venn.diagram(list(pasture_old = pasture_old_AMRs_collapsed$AMR.Gene.Family, pasture_recent = pasture_recent_AMRs_collapsed$AMR.Gene.Family, 
                       urban = urban_AMRs_collapsed$AMR.Gene.Family, natural = natural_AMRs_collapsed$AMR.Gene.Family),
                  fill = c("DarkGreen", "olivedrab4","skyblue3", 'darkorange4'),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)

# genes
v <- venn.diagram(list(pasture_old = unique(pasture_old$Best_Hit_ARO), pasture_recent = unique(pasture_recent$Best_Hit_ARO), 
                       urban = unique(urban$Best_Hit_ARO), natural = unique(natural$Best_Hit_ARO)),
                  fill = c("DarkGreen", "olivedrab4","skyblue3", 'darkorange4'),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)

# resistance mechanism
v <- venn.diagram(list(pasture_old = unique(pasture_old$Resistance.Mechanism), pasture_recent = unique(pasture_recent$Resistance.Mechanism), 
                       urban = unique(urban$Resistance.Mechanism), natural = unique(natural$Resistance.Mechanism)),
                  fill = c("DarkGreen", "olivedrab4","skyblue3", 'darkorange4'),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)




# barplot of the different mechanisms per land use

allparts$land.use <- ifelse(allparts$sample %in% c("RC-T1-3","RC-T1-6","RC-T1-9","SC-T1-3","SC-T1-6","SC-T1-9"), "natural",NA)

allparts$sample <- as.factor(allparts$sample)

#add a column of land use to group the barplot
## This should work and it is not for some reason! 
allparts <- 
  mutate(land.use = case_when(
    sample %in% c("11B-T1-3","11B-T1-6","11B-T1-9",
                  "8-T1-3","8-T1-6","8-T1-9") ~"pasture_recent", 
    sample %in% c("Ex45-T1-3","Ex45-T1-6","Ex45-T1-9",
                  "UAB-T1-3","UAB-T1-6","UAB-T1-9") ~"pasture_old",
    sample %in% c("HP-T1-3","HP-T1-6","HP-T1-9",
                  "RP-T1-3","RP-T1-6","RP-T1-9") ~"urban",
    sample %in% c("RC-T1-3","RC-T1-6","RC-T1-9",
                  "SC-T1-3","SC-T1-6","SC-T1-9") ~"natural"), 
    land.use = factor(land.use, levels = c('pasture_recent', 'pasture_old', 'urban', 'natural')))

    
  
ggplot(allparts, aes(x=Resistance.Mechanism, group = sample, color=sample)) +
  geom_bar()+
  theme(legend.position = "none", text = element_text(size=16))







