library(edgeR)
library(limma)
library(clusterProfiler)
library(KEGGREST)
library(pheatmap)
library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(stringr)
library(goseq)
library(biomaRt)



files <- c("data/counts_SRR8695401.txt",
           "data/counts_SRR8695402.txt",
           "data/counts_SRR8695399.txt",
           "data/counts_SRR8695400.txt",
           "data/counts_SRR8695397.txt",
           "data/counts_SRR8695398.txt")

sample_names <- c("control_1", "control_2", "control_3", "anoxia_1", "anoxia_2", "anoxia_3")


count_list <- lapply(files, function(file) {
  data <- read.delim(file, comment.char = "#", header = TRUE)
  data[, ncol(data)]
})

count_matrix <- do.call(cbind, count_list)
colnames(count_matrix) <- sample_names
rownames(count_matrix) <- read.delim(files[1], comment.char = "#", header = TRUE)$Geneid

dge <- DGEList(counts = count_matrix)
dge$samples$group <- factor(c("control", "control", "control", "anoxia", "anoxia", "anoxia"))
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

group <- dge$samples$group
design <- model.matrix(~group)
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, number=Inf, sort.by="P")

limma::plotMA(fit, coef=2, ylim=c(-5, 5), main="MA-plot")

top_genes <- head(order(res$adj.P.Val), 50)
pheatmap(v$E[top_genes, ], cluster_cols=TRUE, cluster_rows=TRUE, scale="row")

sig_genes <- res[which(res$adj.P.Val < 0.05 & abs(res$logFC) > 1), ]
sig_gene_ids <- rownames(sig_genes)

write.csv(sig_genes, file = "significant_genes.csv")

gff <- import.gff("data/GCF_013100865.1_CAS_Tse_1.0_genomic.gff")
genes <- gff[gff$type == "gene"]
df_genes <- as.data.frame(genes)
head(df_genes[, c("seqnames", "start", "end", "strand", "ID", "Name", "gene", "Note")])


res$gene_id <- rownames(res)
df_genes$gene <- as.character(df_genes$gene)
res$gene_id <- as.character(res$gene_id)
res_annotated <- left_join(res, df_genes, by = c("gene_id" = "gene"))
res_annotated_clean <- res_annotated %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = "; "))))
write.csv(res_annotated_clean, file = "DEG_with_annotations.csv", row.names = FALSE)
head(res_annotated[, c("gene_id", "logFC", "adj.P.Val", "Name", "Note", "seqnames", "start", "end")])

res$threshold <- as.factor(abs(res$logFC) > 1 & res$adj.P.Val < 0.05)

ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value")

pca <- prcomp(t(v$E), scale. = TRUE)

pca_df <- data.frame(Sample = rownames(pca$x),
                     PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     Group = group)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE) +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2") +
  theme_minimal()

plotMDS(dge, col = as.numeric(group), main = "MDS plot")
legend("topright", legend = levels(group), col = 1:2, pch = 20)

top_genes <- head(order(res$adj.P.Val), 50)

pheatmap(v$E[top_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Top 50 DEG Heatmap")

colnames(df_genes)
head(df_genes$Dbxref, 10)