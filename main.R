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

BiocManager::install("genbankr")
library(genbankr)

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

