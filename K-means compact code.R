library(wesanderson)
library(gplots)
library(ComplexHeatmap)

setwd("C:/Users/nateh/OneDrive/Documents/Mei Lab/Categorizing Cardiac Models/Data Analysis/")


k <- 3

kmeans_result <- kmeans(scaledata, centers = k)
cluster <- kmeans_result$cluster
annotation_df <- data.frame(cluster)

#row_annotation <- 


#Appending cluster data to our data table of interest
t <- cbind(scaledata,Cluster)

write.table(annotation_df, "./OrganoidOnly_genelist.txt")
genelist <- read.table("OrganoidOnly_genelist.txt", sep = "\t", header = T)

n_clusters <- length(unique(t$cluster))
annotation_colors <- c("Cluster 1" = "red", "Cluster 2" = "blue", "Cluster 3" = "purple", "Cluster 4" = "orange")




pheatmap(scaledata,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100), # Color scale
         cluster_rows = T,
         cluster_cols = F,
         annotation_row = annotation_df,
         #annotation_colors = annotation_colors,
         main = "K-Means Clustering of Organoid Counts",
         fontsize_row = 8,
         fontsize_col = 8,
         gaps_col = c(9,12,18),
         legend_labels = "Z-Score"
         )

#Attempt to use ComplexHeatmap to solve annotation issue

tscaledata <- t(scaledata)

columnn_anno <- HeatmapAnnotation(df = genelist)

Heatmap(tscaledata, bottom_annotation = columnn_anno)


########-------- YouTube Example for Complex Heatmaps ---------#########

set.seed(123)
mat = matrix(rnorm(100),10)
rownames(mat) = paste0("R", 1:10)
colnames(mat) = paste0("C", 1:10)
column_ha <- HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
row_ha <- rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))

Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)

















