# Visualization and exploratory analyses of the funtional trait data for the Land Use project
# Meri Touceda-Suarez
# March 2021


# --------------------------------------------------------------------------------------------- 
## SETUP ## 
# --------------------------------------------------------------------------------------------- 

# Load packages
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(metagenomeSeq) #bioconductor
library(gRodon)
library(Biostrings)
library(ggrepel)
library(lme4)


# Set working directory 
setwd("/Volumes/MeriTSHD/landuse_functionaltraits/")

# Loa data: 
ko_table <- read.table("ko_table.txt", header = T, row.names = 1)
ko_description <- read.table("ko_description.txt",  header = T, sep="\t")
ko_level1 <- read.table("ko_table_level1.txt", header = T, sep="\t", row.names = 1)
ko_level2 <- read.table("ko_table_level2.txt", header = T, sep="\t", row.names = 1)
ko_level3 <- read.table("ko_table_level3.txt", header = T, sep="\t", row.names = 1)
unknown <- read.table("ko_unknown.txt", header = T)
ko2level <- read.table("ko2level_Jan2021.txt", header = T, sep="\t" )
metadata <- read.csv("LU_metadata.csv", row.names = 1)




#first let's get rid of the X at the beginning of the sample names
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

# stablish rownames of unknown table as sample names
rownames(unknown) <- colnames(ko_table)


#sort all based on metadat rows
unknown <- unknown[rownames(metadata),]
ko_table <- ko_table[, rownames(metadata)]
ko_level1 <- ko_level1[, rownames(metadata)]
ko_level2 <- ko_level2[, rownames(metadata)]
ko_level3 <- ko_level3[, rownames(metadata)]

#check!!! 
all(rownames(metadata) == colnames(ko_table))
all(rownames(metadata) == colnames(ko_level1))
all(rownames(metadata) == colnames(ko_level2))
all(rownames(metadata) == colnames(ko_level3))
all(rownames(metadata) == rownames(unknown))

#record total KO counts per sample
metadata$total_abundance <- colSums(ko_table)




# --------------------------------------------------------------------------------------------------------------------------------
## FUNCTIONAL RICHNESS 
# --------------------------------------------------------------------------------------------------------------------------------

# Here we are going to plot the functional richness (total number of different KOs) of the samples
#BOX richness, number of different KOs
# calculate richness using package vegan
metadata$rich.known <- specnumber(t(ko_table))

# ------------------ Box plot
ggplot(metadata, aes(y = rich.known, x = type, fill = type))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) + 
  ylab("Number of Different KOs")+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_bw()+
  theme(legend.position = "none")
# low outliers are Ex45 in pasture



# ------------------  Eerror line plot

# 1. calculate summary stats
library(dplyr)
rich.summary <- metadata %>%
  group_by(type) %>%
  summarise(
    sd = sd(rich.known, na.rm = TRUE),
    rich = mean(rich.known)
  )
rich.summary

# 2. Plot
ggplot(rich.summary, aes(x = type, y = rich, ymin = rich-sd, ymax = rich+sd, color = type)) +
  geom_errorbar(width = 0.2) +
  geom_point(size = 4) +
  xlab(NULL) + 
  ylab("Number of Different KOs")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_bw()+
  theme(legend.position = "none")

#same thing than "specnumber"
#ko_pa <- decostand(ko_table, method = "pa")
#metadata$parich <- colSums(ko_pa)

# Model
anova(lm(rich.known~type, data = metadata))



# --------------------------------------------------------------------------------------------------------------------------------
## FUNCTIONAL 'COMMUNITY' COMPOSITION ###
# --------------------------------------------------------------------------------------------------------------------------------


# -------------------------------- Normalization? 
hist(rowSums(ko_table))

genes <- colSums(ko_table)
genes
summary(genes)

# Normalize
### No hacer normalizacion porque el numero de genes es muy parecido
## no esta claro si hace falta normalizar en funcitonal
# YJ: # Normalization (metagenomeSeq)
#MR <- newMRexperiment(t(ko_table))
#p <- cumNormStat(MR)
#MR <- cumNorm(MR, p = p)
#ko_norm <- t(MRcounts(MR, norm = T, log = F))


# -------------------------------- Functional Dissimilarity of KOs

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



# --------------------------------------------------------------------------------------------------------------------------------
## FUNCTIONAL TRAITS ##
# --------------------------------------------------------------------------------------------------------------------------------

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



# --------------------------------------------------------------------------------------------------------------------------------
## UNKNOWNS ##
# --------------------------------------------------------------------------------------------------------------------------------

#record proportion of unknown counts out of total counts per sample
metadata$prop.unknown<- unknown$unannotated_number/unknown$gene_number
metadata$prop.annotated.unknown<- unknown$annotated_unknown_number/unknown$gene_number
metadata$prop.annotated.general<- unknown$annotated_general_count/unknown$gene_number


##--------------------------------------------- Proportion genes
# Boxplot Proportion of unnanotated genes
ggplot(metadata, aes(y = prop.unknown,x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) +
  ylab("Proportion of unnaanotated genes (%)")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()+
  theme(legend.position = "none")

# Boxplot proportion of annotated unknown genes
ggplot(metadata, aes(y = prop.annotated.unknown, x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) +
  ylab("Proportion of genes annotated  as 'unknown' (%) ")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()+ 
  theme(legend.position = "none")

# Boxplot proportion of annotated 'unknown'general' genes
ggplot(metadata, aes(y = prop.annotated.general, x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) +
  ylab("Proportion of genes annotated  as 'general' (%) ")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()+ 
  theme(legend.position = "none")


##--------------------------------------------- Gene counts
# Boxplot of unnanotated counts (elative abundance)
ggplot(metadata, aes(y = unknown$unannotated_count, x = type, fill = type))+ #label = metadata$sample
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  #geom_label_repel()+
  xlab(NULL) +
  ylab("Unnanotated gene counts")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()+ 
  theme(legend.position = "none")

# Boxplot of annotated unknown counts
ggplot(metadata, aes(y = unknown$annotated_unknown_count, x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) +
  ylab("Annotated 'unknown' gene counts")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()+ 
  theme(legend.position = "none")

# Boxplot of annotated ''general' counts
ggplot(metadata, aes(y = unknown$annotated_general_count, x = type, fill = type))+
  geom_boxplot(alpha=0.8)+
  geom_jitter(aes(color = type))+ 
  xlab(NULL) +
  ylab("Annotated 'general' gene counts")+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  theme_classic()


# --------------------------------------------------------------------------------------------------------------------------------
## COMPARE FUNCTIONAL PROFILES ##
# --------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------- Kruskal-Wallis
## Non parametric test, used because of low number of samples
# -------------------------------------------------------------- level 3
# Create a relative abundance table of the ko level 3 table
ko_level3_rel <- decostand(ko_level3, method = 'total', margin = 2)

# run kruskal wallis
kruskal.results3<-(apply(t(ko_level3_rel), 2, function(x){kruskal.test(x, metadata$type)}))
kruskal.results3<- data.frame(do.call(rbind, kruskal.results3))
kruskal.adjusted3 <- as.data.frame(p.adjust(kruskal.results3$p.value, method = 'BH',n = length(kruskal.results3$p.value)))
colnames(kruskal.adjusted3) <- 'pvalue'
kruskal.sig3 <- kruskal.adjusted3[which(kruskal.adjusted3$pvalue <= 0.05),]


# -------------------------------------------------------------- level 2
# Create a relative abundance table of the ko level 2 table
ko_level2_rel <- decostand(ko_level2, method = 'total', margin = 2)

# run kruskal wallis
kruskal.results2<-(apply(t(ko_level2_rel), 2, function(x){kruskal.test(x, metadata$type)}))
kruskal.results2<- data.frame(do.call(rbind, kruskal.results2))
kruskal.adjusted2 <- as.data.frame(p.adjust(kruskal.results2$p.value, method = 'fdr',n = length(kruskal.results2$p.value)))
colnames(kruskal.adjusted2) <- 'pvalue'
kruskal.sig2 <- kruskal.adjusted2[which(kruskal.adjusted2$pvalue <= 0.05),]



# -------------------------------------------------------------------------------------------- DeSeqe
# Create the DeSeq object
dds <- DESeqDataSetFromMatrix(countData = ko_level3,  
                              colData = metadata,           
                              design= ~ type)

# Run the algorithm
diagdds = DESeq(dds, test="Wald", fitType="parametric")

resultsNames(diagdds) ###this are all the conditions it tested

# See summary results and some analytics of the process
res = results(diagdds,  name=c("type_pasture.old_vs_natural"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
mcols(res, use.names=TRUE) #explain what are the columns are
head(res)
summary(res) # only one downregulated, one filtered with low counts
DESeq2::plotMA(res) #scatter plot of log2 fold changes
ggplot(as(res, "data.frame"), aes(x = pvalue)) +   #plot a p-value plot
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

# shrunken thing
resLFC <- lfcShrink(diagdds, coef="type_pasture.old_vs_natural", type="apeglm")
DESeq2::plotMA(resLFC) 

# filter only significant 
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
#sigtab = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
head(sigtab)

###order the table by abundance
sigtab = sigtab[order(-sigtab$baseMean),]
####select only the 50 most abundant ones, otherwise we cannot see anything
dim(sigtab)
sigtab50 = sigtab[1:49,]
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
sigtab50.df <- as.data.frame(sigtab50)


# Plot! 

# Order the output of deseq to show level 2 KOs together, better visuals

## aggregate the reference table to have only unique categories
### (the ko2level table has a row for every KO, but categories have multiple KOs, which makes 
### this table very repetitive, and makes it hard to assign levels based on superior or inferior levels (for example, 
### what is the level 1 of the level 2 category 'Translation))
# (only run once) ko2level_unique <- aggregate(KO~level1+level2+level3, ko2level, FUN = length ) 

# create column in sigdata that stores the level3 dfunctions: 
sigtab50.df$level3 <- rownames(sigtab50.df)

# check that all level 3s are in ko2table
sigtab50.df$level3 %in% as.character(ko2level_unique$level3)
# if it is not the case, remove those that have no match in ko2level
#sigtab50.df<- sigtab50.df[which(sigtab50.df$level3 %in% as.character(ko2level_unique$level3)),]


### loop to assign kolevel1 and 2 to sigtab50 based on the level 3 of each row
sigtab50.df$level2 <- ko2level_unique$level2[which(sigtab50.df$level3 %in% as.character(ko2level_unique$level3))] 
sigtab50.df$level1 <- ko2level_unique$level1[which(sigtab50.df$level3 %in% as.character(ko2level_unique$level3))] 



# sort sigtable50 based on level1 
sigtab50.df.sorted <- sigtab50.df[order(sigtab50.df$level1),]

# make level an ordered factor (so ggplot won't put it back into alphabetical order)
sigtab50.df.sorted$level3<- factor(sigtab50.df.sorted$level3, levels=sigtab50.df.sorted$level3)


ggplot(sigtab50.df.sorted, aes(x = log2FoldChange, y = sigtab50.df.sorted$level3, fill = color)) + 
  geom_col(border = NA) + 
  geom_vline(xintercept =  0, size = 0.3)+
  scale_fill_manual(values=c('darkorange4','olivedrab4'))+
  labs(title=NULL)+
  ylab(NULL)+
  xlab(expression(atop("Difference in average relative abundance", paste("(recent pastures- natural soil, %)"))))+
  #xlab('Difference in average relative abundance (recent pasture - natural soil, %)', )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust = TRUE))+
  theme(legend.position = "none")+
  coord_flip()
#geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5) + 
#geom_errorbar(aes(ymin=Mean - my_sd, ymax=Mean + my_sd))






# --------------------------------------------------------------------------------------------------------------------------------
## PLOT FUNCTIONAL TRAITS ##
# --------------------------------------------------------------------------------------------------------------------------------

# create a dataframe for this
funct_plot = as.data.frame(t(ko_level3))
funct_plot$sample <- row.names(funct_plot)


# --------------------------------------------- Stress:

# 'Mismatch repair' 
# boxplot
ggplot(data = funct_plot, aes(y = `Mismatch repair`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) + 
  ylab('Relative abundance of genes involved in Mismatch Repair')+
  theme_classic()+
  theme(legend.position = "none")


# 'Transcription factors'  (sigma factors and chaperones)
# boxplot
ggplot(data = funct_plot, aes(y = `Transcription factors`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) + 
  ylab('Relative abundance of genes related with Transcription Factors')+
  theme_classic()+
  theme(legend.position = "none")



# --------------------------------------------- Competition

# 'Antimicrobial resisteance'
ggplot(data = funct_plot, aes(y = `Antimicrobial resistance genes`/genes, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) + 
  ylab("Relative abundance of Antimicrobial Resistence Genes")+
  theme_classic()+
  theme(legend.position = "none")



# --------------------------------------------- Nitrogen Metabolism
# Nitrogen Metabolism
ggplot(data = funct_plot, aes(y = `Nitrogen metabolism`/metadata$total_abundance, x = metadata$type, fill = metadata$type))+
  geom_boxplot(alpha = 0.8)+
  geom_jitter(aes(color =  metadata$type))+
  scale_fill_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  scale_color_manual(values=c('darkorange4', "DarkGreen", "olivedrab4", "skyblue3"))+
  xlab(NULL) + 
  ylab("Relative abundance of genes involved in Nitrogen Metabolism")+
  theme_classic()+
  theme(legend.position = "none")





