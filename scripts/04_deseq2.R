## DESeq2: genes at peak bleaching (Oct 2024) associated with later survival, accounting for depth
## Inputs (same folder): acer_counts.tsv (from 03_merge_counts.py), sample_table.csv
## One-time install:
##   install.packages(c("ggplot2", "pheatmap", "BiocManager"))
##   BiocManager::install(c("DESeq2", "apeglm"))
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap) })

counts_file <- "acer_counts.tsv"
meta_file   <- "sample_table.csv"
out <- "deseq2_out"; dir.create(out, showWarnings = FALSE)

## ---- 1. Load and align counts to the sample table ----
cts  <- read.delim(counts_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
stopifnot(all(meta$sample %in% colnames(cts)))          # every kept sample has counts
rownames(meta) <- meta$sample
cts <- as.matrix(cts[, meta$sample]); storage.mode(cts) <- "integer"   # drops excluded samples, matches order
stopifnot(identical(colnames(cts), rownames(meta)))

meta$fate     <- factor(meta$fate, levels = c("Dead", "Alive"))        # reference = Dead
meta$depth    <- factor(paste0(meta$depth_ft, "ft"), levels = c("15ft", "85ft", "100ft", "120ft"))
meta$genotype <- factor(meta$genotype)
cat("Samples:", ncol(cts), "\n"); print(table(depth = meta$depth, fate = meta$fate))

## ---- 2. Fit a model ----
fit <- function(design) {
  dds  <- DESeqDataSetFromMatrix(cts, meta, design = design)
  keep <- rowSums(counts(dds) >= 10) >= 5                # expressed in at least 5 samples
  cat(sum(keep), "of", nrow(dds), "genes kept\n")
  DESeq(dds[keep, ])
}
get_res <- function(dds) {
  res  <- results(dds, name = "fate_Alive_vs_Dead", alpha = 0.05)   # positive LFC = higher in survivors
  resS <- lfcShrink(dds, coef = "fate_Alive_vs_Dead", res = res, type = "apeglm")
  tbl  <- as.data.frame(res); tbl$log2FC_shrunk <- resS$log2FoldChange
  tbl$gene <- rownames(tbl); tbl[order(tbl$padj), ]
}

## Model A: the proposal's design
ddsA <- fit(~ depth + fate)
tabA <- get_res(ddsA)
write.csv(tabA, file.path(out, "modelA_depth_fate_results.csv"), row.names = FALSE)
sigA <- subset(tabA, !is.na(padj) & padj < 0.05)
cat("Model A (~ depth + fate): ", nrow(sigA), "genes at BH padj < 0.05\n")

## Model B: also adjusts for genotype (fragments of one genotype are not independent)
ddsB <- tryCatch(fit(~ genotype + depth + fate), error = function(e) { message("Model B failed: ", e$message); NULL })
if (!is.null(ddsB)) {
  tabB <- get_res(ddsB)
  write.csv(tabB, file.path(out, "modelB_genotype_depth_fate_results.csv"), row.names = FALSE)
  sigB <- subset(tabB, !is.na(padj) & padj < 0.05)
  cat("Model B (~ genotype + depth + fate): ", nrow(sigB), "genes; overlap with A:",
      length(intersect(sigA$gene, sigB$gene)), "\n")
}

## ---- 3. QC: PCA on variance-stabilised counts ----
vsd <- vst(ddsA, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("fate", "depth", "genotype"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p1 <- ggplot(pca, aes(PC1, PC2, colour = fate, shape = depth)) + geom_point(size = 3) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) + theme_bw()
ggsave(file.path(out, "pca_fate_depth.png"), p1, width = 6, height = 5, dpi = 300)
p2 <- ggplot(pca, aes(PC1, PC2, colour = genotype, shape = fate)) + geom_point(size = 3) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) + theme_bw()
ggsave(file.path(out, "pca_genotype.png"), p2, width = 6.5, height = 5, dpi = 300)

## ---- 4. Volcano plot (Model A) ----
v <- tabA[!is.na(tabA$padj), ]
v$status <- ifelse(v$padj < 0.05, ifelse(v$log2FC_shrunk > 0, "Higher in survivors", "Higher in dead"), "n.s.")
p3 <- ggplot(v, aes(log2FC_shrunk, -log10(pvalue), colour = status)) + geom_point(size = 1, alpha = 0.6) +
  scale_colour_manual(values = c("Higher in survivors" = "#1b7837", "Higher in dead" = "#b2182b", "n.s." = "grey70")) +
  labs(x = "log2 fold change, survived vs died (shrunk)", y = "-log10 p-value", colour = NULL) + theme_bw()
ggsave(file.path(out, "volcano_modelA.png"), p3, width = 6, height = 5, dpi = 300)

## ---- 5. Heat map of top significant genes (Model A) ----
if (nrow(sigA) >= 2) {
  top <- head(sigA$gene, 50)
  mat <- assay(vsd)[top, ]; mat <- mat - rowMeans(mat)
  ann <- data.frame(fate = meta$fate, depth = meta$depth, row.names = rownames(meta))
  pheatmap(mat, annotation_col = ann, show_rownames = nrow(mat) <= 50, fontsize_row = 6,
           filename = file.path(out, "heatmap_top_genes_modelA.png"), width = 8, height = 8)
} else message("Fewer than 2 significant genes in Model A; no heat map made.")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))   # versions for your report
