# ARGs analyses, land use project
# Mery Touceda-Suárez
# September 2022


# load libraries
library(tidyverse)
library(ggplot2)
library(stringr)

# load and process data
part1 <- read.table("Downloads/card_part1.out.txt", sep ="\t", header = T, fill = T)
part2 <- read.table("Downloads/card_part2.out.txt", sep ="\t", header = T, fill = T)
part3 <- read.table("Downloads/card_part3.out.txt", sep ="\t", header = T, fill = T)
metadata <- read.table("Downloads/metadata-landuse.txt", sep ="\t", header = T)

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
ggplot(metadata, aes(x = land.use, y = args_number))+
  geom_boxplot()


# next steps: 
# QUESTION : should i do the unique value of each sample, or should I put them all together by land use and then do the unique value?
# 1. match the ARG to the "abundance" of such contig or ORF, to see the most abundant ARGs in each sample
# 2. take the list of ARGs in the different land use types and search their potential origin
unique(allparts[allparts$sample=="11B-T1-3",]$Best_Hit_ARO)

pasture_recent <- allparts[allparts$sample %in% c("11B-T1-3","11B-T1-6","11B-T1-9",
                                                  "8-T1-3","8-T1-6","8-T1-9"),]
length(unique(pasture_recent$Best_Hit_ARO)) # number of different ARGs
unique(pasture_recent$Drug.Class)
unique(pasture_recent$Resistance.Mechanism)
unique(pasture_recent$AMR.Gene.Family)

pasture_old <- allparts[allparts$sample %in% c("Ex45-T1-3","Ex45-T1-6","Ex45-T1-9",
                                                  "UAB-T1-3","UAB-T1-6","UAB-T1-9"),]
length(unique(pasture_old$Best_Hit_ARO))
unique(pasture_old$Drug.Class)
unique(pasture_old$Resistance.Mechanism)
unique(pasture_old$AMR.Gene.Family)

urban <- allparts[allparts$sample %in% c("HP-T1-3","HP-T1-6","HP-T1-9",
                                                  "RP-T1-3","RP-T1-6","RP-T1-9"),]
length(unique(urban$Best_Hit_ARO))
unique(urban$Drug.Class)#different drug classes
length(unique(urban$Drug.Class))#number of different drug classes
unique(urban$Resistance.Mechanism)
unique(urban$AMR.Gene.Family)#different families 
length(unique(urban$AMR.Gene.Family))#number of different families 

natural <- allparts[allparts$sample %in% c("RC-T1-3","RC-T1-6","RC-T1-9",
                                                  "SC-T1-3","SC-T1-6","SC-T1-9"),]
length(unique(natural$Best_Hit_ARO))
unique(natural$Drug.Class)#different drug classes
length(unique(natural$Drug.Class))#number of different drug classes
unique(natural$Resistance.Mechanism)
unique(natural$AMR.Gene.Family)#different families 
length(unique(natural$AMR.Gene.Family))#number of different families
