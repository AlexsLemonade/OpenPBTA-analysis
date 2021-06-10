suppressPackageStartupMessages(library(tidyverse))

#------------ Create output directories --------------------
table_outdir <- file.path('results', paste0('deseq2_rle', '_normalized'))
dir.create(table_outdir, showWarnings = FALSE)

plot_outdir <- file.path('plots', paste0('deseq2_rle', '_normalized'))
dir.create(plot_outdir, showWarnings = FALSE)
#------------ Read and pre-process data --------------------
print('Read count matrices...')
cnt_mat <- readRDS('../../data/gene-counts-rsem-expected_count-collapsed.rds')

# 6 samples have both poly-A and stranded RNA-seq libraries
# Index | Kids_First_Biospecimen_ID | sample_id | experimental_strategy | 
#   RNA_library | cohort
# -- | -- | -- | -- | -- | --
# 1 | BS_HE0WJRW6 | 7316-1455 | RNA-Seq | stranded | CBTN
# 2 | BS_HWGWYCY7 | 7316-1455 | RNA-Seq | poly-A | CBTN 
# 3 | BS_SHJA4MR0 | 7316-161 | RNA-Seq | stranded | CBTN 
# 4 | BS_X0XXN9BK | 7316-161 | RNA-Seq | poly-A | CBTN 
# 5 | BS_FN07P04C | 7316-255 | RNA-Seq | stranded | CBTN 
# 6 | BS_W4H1D4Y6 | 7316-255 | RNA-Seq | poly-A | CBTN
# 7 | BS_8QB4S4VA | 7316-536 | RNA-Seq | stranded | CBTN
# 8 | BS_QKT3TJVK | 7316-536 | RNA-Seq | poly-A | CBTN
# 9 | BS_7WM3MNZ0 | A16915 | RNA-Seq | poly-A | PNOC003
# 10 | BS_KABQQA0T | A16915 | RNA-Seq | stranded | PNOC003
# 11 | BS_68KX6A42 | A18777 | RNA-Seq | poly-A | PNOC003
# 12 | BS_D7XRFE0R | A18777 | RNA-Seq | stranded | PNOC003
counts <- cbind(
  cnt_mat[, c("BS_HE0WJRW6", "BS_SHJA4MR0", "BS_FN07P04C",
              "BS_8QB4S4VA", "BS_KABQQA0T", "BS_D7XRFE0R")],
  cnt_mat[, c("BS_HWGWYCY7", "BS_X0XXN9BK", "BS_W4H1D4Y6",
              "BS_QKT3TJVK", "BS_7WM3MNZ0", "BS_68KX6A42")]
)
group <- factor(c(rep("stranded", 6), rep("polya", 6)))

counts_object <- edgeR::DGEList(counts = counts, group = group)
counts_object <- counts_object[
  edgeR::filterByExpr(counts_object), , keep.lib.sizes=FALSE]
#------------ Run DESeq2 nbinomWaldTest ------------------------
print(paste0('Run differential gene expression DESeq2 nbinomWaldTest on',
             ' poly-A vs stranded RNA-seq...'))
# DESeq2 requires integer matrix as read count matrix
round_cnt_mat <- round(counts_object$counts)
suppressMessages(
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round_cnt_mat,
    colData = data.frame(group),
    design = ~ group)
)
# Read housekeeping genes
hkgenes <- read.csv(file.path('input', 'Housekeeping_GenesHuman.csv'),
                    sep = ";", stringsAsFactors = FALSE)
hkgenes <- unique(hkgenes$Gene.name)
hkgenes <- intersect(rownames(counts_object), hkgenes)
# From DESeq2 documentation:
# DESeq2 performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
dds <- DESeq2::estimateSizeFactors(dds)
dds <- dds[hkgenes, ]
suppressMessages(
  dds <- DESeq2::estimateDispersions(dds)
)
dds <- DESeq2::nbinomWaldTest(dds)
# DESeq2::counts(dds, normalized=TRUE)[1:5, 1:5]
deseq2_res <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
deseq2_res_df <- data.frame(deseq2_res)
deseq2_res_df <- deseq2_res_df[order(deseq2_res_df$pvalue), ]

deseq2_out_df <- deseq2_res_df
stopifnot(identical(colnames(deseq2_out_df),
                    c("baseMean", "log2FoldChange", "lfcSE",
                      "stat", "pvalue", "padj")))
colnames(deseq2_out_df) <- c("mean_counts", "stranded_over_polya_log2FC",
                             "stranded_over_polya_log2FC_standard_error",
                             "statistic", "PValue", "BH FDR")
write.csv(
  deseq2_out_df,
  file=file.path(table_outdir,
                 'stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv'))

# plot and save p-value histogram
# some DESeq2 p-values are NAs due to all-0 normalized read counts
deseq2_plot_df <- deseq2_res_df[, c('pvalue', 'padj')]
deseq2_plot_df[is.na(deseq2_plot_df)] <- 1
p <- ggplot(deseq2_plot_df, aes(x=pvalue)) +
  geom_histogram(binwidth = 0.05, center = 0.025) +
  theme_classic() +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  xlab('stranded vs poly-A RNA-seq DGE RLE nbinomWaldTest p-value') +
  ylab('Gene count') +
  ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                 'differential gene expression RLE normalized\n',
                 'DESeq2 nbinomWaldTest p-values\n',
                 sum(deseq2_plot_df$pvalue < 0.05),
                 ' genes have p-value < 0.05\n',
                 sum(deseq2_plot_df$pvalue >= 0.05),
                 ' genes have p-value >= 0.05\n',
                 sum(deseq2_plot_df$padj < 0.05),
                 ' genes have BH FDR < 0.05\n',
                 sum(deseq2_plot_df$padj >= 0.05),
                 ' genes have BH FDR >= 0.05')) +
  theme(text = element_text(size=15))

ggsave(
  file.path(
    plot_outdir,
    'stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png'),
  dpi = 300, plot = p, width = 8, height = 7)
print('Done.')
