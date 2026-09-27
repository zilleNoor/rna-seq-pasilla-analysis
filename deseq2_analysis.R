library(DESeq2)
library(pasilla)
library(ggplot2)

# Load counts matrix
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
counts <- read.table(pasCts, header=TRUE, row.names=1)

# Load sample metadata
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)
metadata <- read.csv(pasAnno, row.names=1)

# metadata sample names had an extra "fb" suffix that counts didn't have — stripped it so both tables match
rownames(metadata) <- sub("fb$", "", rownames(metadata))
metadata <- metadata[colnames(counts), ]

# Build DESeq2 dataset and run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

# Filter to significant genes (padj < 0.05), sorted by significance
res_sig <- res[which(res$padj < 0.05), ]
res_sig <- res_sig[order(res_sig$padj), ]
write.csv(as.data.frame(res_sig), "significant_genes.csv")

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot: Treated vs Untreated",
       x = "log2 Fold Change",
       y = "-log10(p-value)")

ggsave("volcano_plot.png", width = 6, height = 5, dpi = 300)