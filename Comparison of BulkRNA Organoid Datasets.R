library(ggplot2)
library(DESeq2)
library(org.Hs.eg.db)
library(sva)
library(tidyverse)
library(clusterProfiler)
library(ggthemes)
library(enrichplot)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(rgl)
library(scatterplot3d)

##### Loading datasets into RStudio #####

setwd("C:/Users/nateh/OneDrive/Documents/Mei Lab/Categorizing Cardiac Models/Data Analysis/")

colnames <- read.table("Colnames/colnames_organoids.txt", header = T)
counts <- as.data.frame(read.csv("Counts/counts_organoid_only.csv", header = T))

###------- Changing the GeneID's to Symbols- just makes for clearer data

columns(org.Hs.eg.db)
symbol <- counts$Geneid

test <- mapIds(org.Hs.eg.db, as.character(symbol), 'SYMBOL', 'ENTREZID')
test2 <- as.data.frame(test[1:length(test)])
colnames(test2)[1] <- "Symbol"
Best <- merge(test2, counts, by.x = 0 , by.y = 1 , all.x = TRUE)
Best$Row.names <- NULL
colnames(Best)[1] <- "Symbol"
counts <- Best
write.table(Best, "./counts_Symbol.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###------- Formatting Counts for DESeq Use

counts <- aggregate(.~Symbol, data = counts, max)
counts <- counts[!is.na(counts$Symbol),]
rownames(counts) <- counts$Symbol
counts$Symbol <- NULL


#Making a DESeq object, sorting by group
dds <- DESeqDataSetFromMatrix( countData = counts,
                               colData = colnames,
                               design = ~Group)

#Rlog normalizaiton of the data
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <-rlog(dds, blind = FALSE)

rldmatrix <- as.matrix(assay(rld))
rldmatrix <- as.data.frame(rldmatrix)
rlogquant <- rldmatrix
rownames(rlogquant) <- rownames(rldmatrix)
colnames(rlogquant) <- colnames(rldmatrix)

write.table(rlogquant, "./DESEQ_counts_hCOs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#The logquant is now transposed (not sure why), then scaled using the Z score
transposerlogquant <- t(rlogquant)
tscalerlogquant <- scale(transposerlogquant) #scale using Zscore
zScoreRlog <- t(tscalerlogquant)
zScoreRlog <- zScoreRlog[complete.cases(zScoreRlog),]

rvquant <- rowVars(as.matrix(rlogquant))
select <- order(rvquant, decreasing=TRUE)[seq_len(min(500, length(rvquant)))]

pcaquant <- prcomp(t(rlogquant[select,])) #prcomp does the PCA analysis and returns a matrix of numbers rather than a graph 
df_quant <- as.data.frame(pcaquant$x)
df_quant$Group <- sapply( strsplit(as.character(colnames$Group), "_"), "[[", 1)
df_quant$Subgroup <- colnames$Subgroup


##### ---------- PCA Plot ---------- ######
autoplot(pcaquant,
         data = colnames,
         colour = "Group",
         frame = T) +
    ggtitle("Principal Component Analysis of hCOs") + 
    xlab("PC1: 48% Variance") + 
    ylab("PC2: 41% Variance") + 
    theme_gray(base_size = 14) +
    scale_color_manual(values = c("#E7298A", "#66A61E", "#666666", "#A6761D")) +
    scale_fill_manual(values = c("#E7298A", "#66A61E", "#666666", "#A6761D"))


#Creating a gene-ranking object
Gene <- rownames(pcaquant$rotation)
Gene <- as.data.frame(Gene)
Gene$Rank <- pcaquant$rotation[,2]
write.table(Gene, "PC2rankedgenes_Organoids.rnk", sep = "\t", quote = FALSE, row.names = FALSE)
