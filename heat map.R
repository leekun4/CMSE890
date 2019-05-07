source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(tidyverse)
library(pheatmap)

top_point_five_genes <- read.csv("top_0.5_percent_var_genes.csv")
top_five_genes_flipped <- scale(t(top_point_five_genes), center = T, scale = T)
pheatmap(top_five_genes_flipped)

top_point_ten_genes <- read.csv("top_10_percent_var_genes.csv")
top_ten_genes_flipped <- scale(t(top_point_ten_genes), center = T, scale = T)
pheatmap(top_ten_genes_flipped)

gene_expression <- "FPKM_expression_profile.txt"
my_data <- read.table(gene_expression, sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)

exp_matrix_nonzero_top20percvargenes
dim(exp_matrix_nonzero_top20percvargenes)
class(exp_matrix_nonzero_top20percvargenes)
nrow(exp_matrix_nonzero_top20percvargenes)
ncol(exp_matrix_nonzero_top20percvargenes)
nrow(exp_matrix_nonzero_top10percvargenes)
class(exp_matrix_nonzero_top10percvargenes)
head(exp_matrix_nonzero_top10percvargenes)
nrow(top_point_five_genes)

#exporting a file as a csv 
write.csv(as.data.frame(exp_matrix_nonzero_top0.5percvargenes), 
          file="top_0.5_percent_var_genes.csv")

Heatmap(top_point_five_genes)

#flip rows and columns around
top_five_genes_flipped <- t(top_point_five_genes)
Heatmap(top_five_genes_flipped)

#not clustering
Heatmap(top_point_five_genes, cluster_columns=FALSE)

fontsize <- 0.6

Heatmap(top_point_five_genes,
        cluster_columns=FALSE,
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=fontsize),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        km=2) # number of clusters you want



