install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
library(DESeq2)
library(pasilla)
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
counts <- read.table(pasCts, header=TRUE, row.names=1)
head(counts)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)
metadata <- read.csv(pasAnno, row.names=1)
head(metadata)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition )
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition )
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = metadata,
design = ~condition)
colnames(counts)
rownames(count)
rownames(metadata)
metadata <- metadata[colnames(counts),]
rownames(metadata)
colnames(counts)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = metadata,
design = ~condition)
identical(colnames(counts), rownames(metadata))
colnames(counts)
rownames(metadata)
rownames(metadata) <- sub("fb$","",rownames(metadata))
identical(colnames(counts), rownames(metadata))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = metadata,
design = ~condition)
dds
res<- results(dds)
head(res)
res_sig <- res[which(res$padj < 0.05), ]
res_sig <- res_sig[order(res_sig$padj), ]
head(res_sig)
write.csv(as.data.frame(res_sig), "significant_genes.csv")
head(res)
res_sig <- res[which(res$padj < 0.05), ]
res_sig <- res_sig[order(res_sig$padj), ]
write.csv(as.data.frame(res_sig), "significant_genes.csv")
nrow(res_sig)
library(ggplot2)
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
geom_point(alpha = 0.5) +
theme_minimal() +
labs(title = "Volcano Plot: Treated vs Untreated",
x = "log2 Fold Change",
y = "-log10(p-value)")
ggsave("volcano_plot.png", width = 6, height = 5, dpi = 300)
