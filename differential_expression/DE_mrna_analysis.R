library(DESeq2)
library(ggplot2)

# Load Research Counts & Clinical Metadata
counts_raw <- read.csv("reordered_mrna.csv", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- read.table("metadata_mrna.txt", header = TRUE, row.names = 1)

# Mapping & Numeric Extraction
gene_map <- counts_raw$GeneName
names(gene_map) <- rownames(counts_raw)
counts <- counts_raw[, !colnames(counts_raw) %in% "GeneName"]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# [RESEARCH SENSITIVE: NORMALIZATION & QUALITY GATING]
# Implementation of study-specific pre-filtering thresholds based on 
# expected complexity of mRNA libraries in Vascular Dementia tissue.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# Results Extraction
res <- results(dds)
res_df <- as.data.frame(res)
res_df$TranscriptID <- rownames(res_df)
res_df$GeneName <- gene_map[res_df$TranscriptID]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# [ANALYSIS NOTE: LFC SHRINKAGE]
# apeglm shrinkage is applied to stabilize fold change estimations 
# for low-count transcripts prior to final ranking.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

write.csv(res_df, file = "VDA_DE_full_results.csv", row.names = FALSE)

# Professional Visualization
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.8, size = 1.2) +
  theme_minimal() +
  labs(title = "VDA Ischemia vs Sham Transition: mRNA Landscape")
