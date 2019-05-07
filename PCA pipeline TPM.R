library(tidyverse)

data.matrix <- read_csv("~/Desktop/Materials for research/DESeq2/FPKM_Gene_Expression.csv") #, row.names = 1)
data.matrix

fpkm <- as.matrix(data.matrix[,-1])

tpm <- apply(fpkm, 2,
              function(x) {
                log2(x*(10**6)/sum(x,na.rm = T) + 1)
              })

data.matrix <- cbind(data.matrix[,1],
                     tpm) %>% as_tibble()

data.matrix %>%
  group_by(X1) %>%
  summarize(n()) %>%
  arrange(desc(`n()`))

gene_tpm_sum <- apply(tpm, 1, sum)
tpm_nonzero <- tpm[(gene_tpm_sum > 0),]

gene_tpm_var <- apply(tpm_nonzero, 1, var)
tpm_nonzero_top50percvargenes <- tpm_nonzero[(gene_tpm_var > quantile(gene_tpm_var, probs = 0.5)),]
tpm_nonzero_top20percvargenes <- tpm_nonzero[(gene_tpm_var > quantile(gene_tpm_var, probs = 0.8)),]
tpm_nonzero_top10percvargenes <- tpm_nonzero[(gene_tpm_var > quantile(gene_tpm_var, probs = 0.9)),]
tpm_nonzero_top1percvargenes <- tpm_nonzero[(gene_tpm_var > quantile(gene_tpm_var, probs = 0.99)),]
tpm_nonzero_top0.5percvargenes <- tpm_nonzero[(gene_tpm_var > quantile(gene_tpm_var, probs = 0.995)),]

tpm_nonzero_top0.5percvargenes %>%
  as.data.frame() %>%
  write_csv("top_0.5_percent_var_genes.csv")

tpm_nonzero_top10percvargenes %>%
  as.data.frame() %>%
  write_csv("top_10_percent_var_genes.csv")

pca <- prcomp(t(tpm_nonzero_top0.5percvargenes), scale=TRUE)


pca <- prcomp(t(tpm_nonzero_top50percvargenes), scale=TRUE)
pca <- prcomp(t(exp_matrix_nonzero), scale=TRUE)



# data.working <- scale(data.matrix[c("M3B5","M10B2", "M9B3", "M1B12", "M3B2", "M3B1",
#                                "M3B6", "M10B3", "M9B2", "M1B2", "M1B4", "M1B3")])
# pca <- prcomp(t(data.working), scale=TRUE)
# plot(pca$x[,1], pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x),
                       PC1=pca$x[,1],
                       PC2=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=PC1, y=PC2, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA analysis")

loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes
pca$rotation[top_10_genes,1]
