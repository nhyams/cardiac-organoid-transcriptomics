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

#For this Gene Panel, we load in our logfoldchange data and colnames files *NOTE- this data is pulled AFTER DESeq and Rlog transformations

setwd("C:/Users/nateh/OneDrive/Documents/Mei Lab/Categorizing Cardiac Models/Data Analysis/")
GenePanel <- as.data.frame(read.csv("counts/Gene_Panel_Data.csv", header = T))
GenePanelcolnames <- read.table("Colnames/colnames_panel.txt")

col_names <- GenePanel[,1]
GenePanel <- GenePanel[,-1]
rownames(GenePanel) <- col_names

#Setting up annotation dataframes for the Gene Panel

#Colors assigned to annotations
anno_color <- list('Biological_Process' = c("Calcium Handling" = "#1f77B4", "Atrial CMs" = "darkred", "Ventricular CMs" = "wheat", "CM Maturation" = "#66A61E", "Metabolism" = "#7570B3"), 
                     'sample' = c("2D CM" = "#1B9E77", "3D CM" = "#D95F02", "EB hCO" = "#E7298A", "Primary Stroma hCO" = "#666666", "Extended Culture hCO" = "#66A61E",
                                  "HeartDyno hCO" = "#A6761D", "Fetal Heart" = "#E6AB02", "Adult Heart" = "#7570B3"))

#Biological Processes assigned to certain genes
BP_anno <- data.frame(Biological_Process = rep(c("Calcium Handling", "Atrial CMs", 
                                          "Ventricular CMs", "CM Maturation", "Metabolism"), c(6, 3, 3, 5, 2)))

row.names(BP_anno) <- colnames(GenePanel)

#Groups of samples
sample_group <- data.frame(sample = rep(c("2D CM", "3D CM", 
                                          "EB hCO", "Primary Stroma hCO", "Extended Culture hCO", 
                                          "HeartDyno hCO", "Fetal Heart", "Adult Heart"), c(3, 3, 9, 3, 6, 4, 4, 3)))
row.names(sample_group) <- rownames(GenePanel)


#Generating a Heatmap of Differentially Expressed Genes of Interest
pheatmap(GenePanel,
         cluster_rows = F,
         cluster_cols = F,
         annotation_colors = anno_color,
         annotation_col = BP_anno,
         annotation_row = sample_group,
         gaps_row = c(3, 6, 15, 18, 24, 28, 32),
         gaps_col = c(6, 9, 12, 17),
         show_rownames = F,
         show_colnames = T,
         cellwidth = 15,
         border_color = "grey")
