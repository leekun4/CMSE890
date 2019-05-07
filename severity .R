library(DESeq2)
library(ggplot2)

#count matrix input
cts <- read.csv("SeverityCounts.csv", sep = ",")
coldata <- read.csv("Infection_severity.csv", sep = ",", row.names = 1)
head(cts,2)
coldata

#making sure they are in order
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#constructing DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Conditions)
dds

#factor levels
dds$Conditions <- factor(dds$Conditions, levels = c("High","Intermediate", "Low"))
dds$Conditions <- relevel(dds$Conditions, ref = "Low")

#differential expression analysis
dds <- DESeq(dds)
res <-results(dds)
res

resultsNames(dds)

#exporting DESeq2 results without any filtering
write.csv(as.data.frame(res), 
          file="Severity_genes.csv")

#apeglm method for effect size shrinkage
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm", version = "3.8")
library(apeglm)

#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Conditions_Low_vs_Intermediate_vs_High", type="apeglm")
resLFC

#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]
summary(res)

#adjusted p-values less than 0.05
sum(res$pvalue < 0.05, na.rm=TRUE)

#MA-plot (red dot means if the adjusted p value is les than 0.1)
plotMA(res, ylim=c(-2,2))

#shrunken log 2 fold changes
plotMA(resLFC, ylim=c(-2,2))

#plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")

#plot counts ggplot
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#exporting results to CSV files
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")

colData(dds)

#extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#download vsn
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("vsn", version = "3.8")

#effects of transformations on the variance
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

#heatmap of the count matrix
install.packages("pheatmap")
library("pheatmap")
library("ggplot2")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds) [,c("Condition")])
rownames(df) <-colnames(select)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#PCA analysis
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

plotPCA(vsd, intgroup=c("Condition"))

#likelihood ratio test 
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

#outliers box plot
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

?pch
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-7,5)))
#adding colors to volcano plot
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, padj>.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
