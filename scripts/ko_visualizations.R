# Visualization and exploratory analyses of the funtional trait data for the Land Use project
# Meri Touceda-Suarez
# March 2021


##################################################################################################################################
#################################################### SETUP####################################################################
##################################################################################################################################

# Load packages
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(metagenomeSeq) #bioconductor
library(gRodon)
library(Biostrings)



setwd("/Volumes/MeriTSHD/landuse_functionaltraits/")
# Loa data: 
ko_table <- read.table("ko_table.txt", header = T, row.names = 1)
ko_description <- read.table("ko_description.txt",  header = T, sep="\t")
ko_level1 <- read.table("ko_table_level1.txt", header = T, sep="\t", row.names = 1)
ko_level2 <- read.table("ko_table_level2.txt", header = T, sep="\t", row.names = 1)
ko_level3 <- read.table("ko_table_level3.txt", header = T, sep="\t", row.names = 1)
unknown <- read.table("ko_unknown.txt", header = T)
ko2level <- read.table("ko2level_Jan2021.txt", header = T, sep="\t" )



#first let's get rid of the X at the beginning
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

destroyX(ko_table)
destroyX(ko_level1)
destroyX(ko_level2)
destroyX(ko_level3)


rownames(unknown) <- colnames(ko_table)



#make a quick metadata to know what samples are urban/pasture/natural
metadata <- data.frame(matrix(vector(), 24, 2,
                              dimnames=list(c(), c("sample", "type"))))
metadata$sample <- rownames(unknown)
metadata$type <- c("pasture-recent", "pasture-recent","pasture-recent", "pasture-recent","pasture-recent", "pasture-recent","pasture-old", "pasture-old","pasture-old", "urban", 
                   "urban", "urban", "natural", "natural", "natural", "urban", "urban", "urban", "natural", "natural", "natural", 
                   "pasture-old", "pasture-old", "pasture-old")



## Incluir columna que sea site! 

#sort all based on type
metadata <- metadata[order(metadata$type),]
unknown <- unknown[order(metadata$type),]
ko_table <- ko_table[, order(metadata$type)]
ko_level1 <- ko_level1[, order(metadata$type)]
ko_level2 <- ko_level2[, order(metadata$type)]
ko_level3 <- ko_level3[, order(metadata$type)]

metadata$site <-c("Rose Canyon","Rose Canyon","Rose Canyon","Sabino Canyon","Sabino Canyon", "Sabino Canyon", 
                  "Ex45", "Ex45", "Ex45", "UAB","UAB", "UAB","11", "11", "11", "8", "8", "8", "Himmel Park", 
                  "Himmel Park","Himmel Park", "Reid Park","Reid Park","Reid Park")

#record total KO counts per sample
metadata$total_abundance <- colSums(ko_table)

rownames(metadata) <- metadata$sample




##################################################################################################################################
#################################################### RICHNESS (alpha diversity) ##################################################
##################################################################################################################################

#BOX richness, number of different KOs
metadata$rich <- specnumber(t(ko_table))
ggplot(metadata, aes(y = rich, x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) + 
  ylab("Number of Different Genes")+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_bw()+
  theme(legend.position = "none")
# low outliers are Ex45 in pasture


## hacer un error barplot en vez de box

#same thing than "specnumber"
#ko_pa <- decostand(ko_table, method = "pa")
#metadata$parich <- colSums(ko_pa)


##################################################################################################################################
#################################################### EVENNESS (Beta diversity) ####################################################
##################################################################################################################################

# Normalize
hist(rowSums(ko_table))

genes <- colSums(ko_table)
genes
summary(genes)

### No hacer normalizacion porque el numero de genes es muy parecido
## no esta claro si hace falta normalizar en funcitonal
# YJ: # Normalization (metagenomeSeq)
#MR <- newMRexperiment(t(ko_table))
#p <- cumNormStat(MR)
#MR <- cumNorm(MR, p = p)
#ko_norm <- t(MRcounts(MR, norm = T, log = F))


ko_bray <- vegdist(t(ko_table), method="bray")
ko_nmds <- metaMDS(ko_bray, k=2, try = 100) 
metadata$Axis01 = ko_nmds$points[,1]
metadata$Axis02 = ko_nmds$points[,2]

ko_nmds$stress 

library(ggforce)
library(concaveman)

ggplot(metadata, aes(Axis01, Axis02))+
  geom_point(aes(color=type), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=type, color = type))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"), name = "Environment")+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"), name = "Environment")+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  #scale_alpha_manual(values=c("Grazed"=0.4, "Ungrazed"=0.9))+
  theme_bw(base_size = 15)


# permanova
adonis(ko_bray~metadata$type)
adonis(ko_bray~metadata$site)

adonis(ko_bray~type/site, strata=metadata$type, data= metadata)
1 - 0.28486  # 0.71514 (1-residuals)


# el type explica el 49% y el site otro 22%

# Efecto site 


##################################################################################################################################
#################################################### FUNCTIONAL TRAITS ###########################################################
##################################################################################################################################

metadata <- metadata[order(metadata$sample),]
ko3_prop <- ko_level3 / metadata$total_abundance

ko1_prop <- ko_level1 / metadata$total_abundance

ggplot(ko1_prop)+
  geom_boxplot(aes(color=metadata$type))

## Do the proportion and then plot KOlevel3 grouped nby KOlevel1


# Richness and Ordination ()
#Ordination: normalize, do bray curtis and then nmds



# For KO levels: 
# organize level 3 levels into categories (using KO2level table)



#DeSeq2 (negative binomial) para los functional (YJ)
# leerme paper de YJ del ISME



##################################################################################################################################
#################################################### UNKNOWNS ####################################################################
##################################################################################################################################

#record proportion of unknown counts out of total counts per sample
metadata$u.relabund <- unknown[, 6] / metadata[, 3]

#plot that for the Natural/Pasture/Urban sites
ggplot(metadata, aes(y = u.relabund, group =type, color = type))+
  ylab("Proportion of annotated unknown counts")+
  geom_boxplot()+
  theme_classic()




##################################################################################################################################
#################################################### COMPARE FUNCTIONAL PROFILES #################################################
##################################################################################################################################

# Comparison of functional profiles
# Need non parametric tests because low sample number
# Use wilcox.test()
# Kruskal-Wallis Rank Sum test: kruskal.test

##################################################################################################################################
#################################################### PLOT FUNCTIONAL TRAITS ######################################################
##################################################################################################################################

# create a dataframe for this
funct_plot = as.data.frame(t(ko_level3))
funct_plot$sample <- row.names(funct_plot)



################## STRESS: 

# 'Mismatch repair' 
# boxplot
ggplot(data = funct_plot, aes(y = `Mismatch repair`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()



# 'Transcription factors'  (sigma factors and chaperones)
# boxplot
ggplot(data = funct_plot, aes(y = `Transcription factors`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()



################## Competition

# Antimicrobial resisteance
ggplot(data = funct_plot, aes(y = `Antimicrobial resistance genes`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) + 
  ylab("Antimicrobial Resistence Genes")+
  theme_bw()+
  theme(legend.position = "none")





################## Growth






