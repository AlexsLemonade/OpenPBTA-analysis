suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(tidyverse))


output_deseq2_res_df <- function(res_df, save_colnames, path) {
    out_df <- res_df
    stopifnot(identical(colnames(out_df),
                        c("baseMean", "log2FoldChange", "lfcSE",
                          "stat", "pvalue", "padj")))
    colnames(out_df) <- save_colnames
    write.csv(out_df, file=path)
}


deseq2_pvals_histogram <- function(res_df, xlab, ylab, title) {
    stopifnot(all(c('pvalue', 'padj') %in% colnames(res_df)))
    p <- ggplot(res_df, aes(x=pvalue)) +
        geom_histogram(binwidth = 0.05, center = 0.025) +
        theme_classic() +
        scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(paste0(title, '\n',
                       sum(res_df$pvalue < 0.05),
                       ' genes have p-value < 0.05\n',
                       sum(res_df$pvalue >= 0.05),
                       ' genes have p-value >= 0.05\n',
                       sum(res_df$padj < 0.05),
                       ' genes have BH FDR < 0.05\n',
                       sum(res_df$padj >= 0.05),
                       ' genes have BH FDR >= 0.05')) +
        theme(text = element_text(size=15))
    return(p)
}


#------------ Create output directories ----------------------------------------
table_outdir <- 'results'
dir.create(table_outdir, showWarnings = FALSE)

plot_outdir <- 'plots'
dir.create(plot_outdir, showWarnings = FALSE)


#------------ Read and pre-process data ----------------------------------------
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


#------------ Run DESeq2 nbinomWaldTest ----------------------------------------
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
# From DESeq2 documentation:
# DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
suppressMessages(
    dds <- DESeq2::DESeq(dds)
)

# Export results and save as CSV
deseq2_res <- DESeq2::results(dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
deseq2_res_df <- data.frame(deseq2_res)
deseq2_res_df <- deseq2_res_df[order(deseq2_res_df$pvalue), ]

output_deseq2_res_df(
    deseq2_res_df,
    save_colnames = c("mean_counts", "stranded_over_polya_log2FC",
                      "stranded_over_polya_log2FC_standard_error",
                      "statistic", "PValue", "BH FDR"),
    path = file.path(table_outdir,
                     'stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv')
)

# plot and save p-value histogram
# some DESeq2 p-values are NAs due to all-0 normalized read counts
deseq2_plot_df <- deseq2_res_df[, c('pvalue', 'padj')]
deseq2_plot_df[is.na(deseq2_plot_df)] <- 1
p <- deseq2_pvals_histogram(
    deseq2_plot_df,
    'stranded vs poly-A RNA-seq DGE RLE nbinomWaldTest p-value',
    'Gene count',
    paste0('Histogram of stranded vs poly-A RNA-seq\n',
           'differential gene expression RLE normalized\n',
           'DESeq2 nbinomWaldTest p-values'))


ggsave(
    file.path(
        plot_outdir,
        'stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png'),
    dpi = 300, plot = p, width = 8, height = 7)


#------------ Run RUVSeq batch correction and DESeq2 nbinomWaldTest ------------
print(paste0('Run differential gene expression DESeq2 nbinomWaldTest ',
             'with RUVSeq estimated batch effect',
             'on poly-A vs stranded RNA-seq...'))
seq_expr_set <- newSeqExpressionSet(
    round_cnt_mat,
    phenoData = data.frame(group, row.names = colnames(round_cnt_mat)))

seg_df <- read.csv(
    'input/uqpgq2_normalized_stranded_vs_polya_stably_exp_genes.csv',
    stringsAsFactors = FALSE, row.names = 1)

ruvg_res <- RUVg(seq_expr_set, seg_df$gene, k=1)

suppressMessages(
    ruv_dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts(ruvg_res),
        colData = pData(ruvg_res),
        design = ~ W_1 + group)
)
# From DESeq2 documentation:
# DESeq2::DESeq performs a default analysis through the steps:
# - estimation of size factors: estimateSizeFactors
# - estimation of dispersion: estimateDispersions
# - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
suppressMessages(
    ruv_dds <- DESeq2::DESeq(ruv_dds)
)

# Export results and save as CSV
ruv_deseq2_res <- DESeq2::results(
    ruv_dds, cooksCutoff = FALSE, pAdjustMethod = 'BH')
ruv_deseq2_res_df <- data.frame(ruv_deseq2_res)
ruv_deseq2_res_df <- ruv_deseq2_res_df[order(ruv_deseq2_res_df$pvalue), ]

output_deseq2_res_df(
    ruv_deseq2_res_df,
    save_colnames = c("mean_counts", "stranded_over_polya_log2FC",
                      "stranded_over_polya_log2FC_standard_error",
                      "statistic", "PValue", "BH FDR"),
    path = file.path(
        table_outdir,
        'stranded_vs_polya_dge_ruvg_k1_deseq2_nbinom_wald_test_res.csv')
)



# plot and save p-value histogram
# some DESeq2 p-values are NAs due to all-0 normalized read counts
ruv_deseq2_plot_df <- ruv_deseq2_res_df[, c('pvalue', 'padj')]
ruv_deseq2_plot_df[is.na(ruv_deseq2_plot_df)] <- 1
p <- deseq2_pvals_histogram(
    ruv_deseq2_plot_df,
    'stranded vs poly-A RNA-seq DGE RLE RUVg nbinomWaldTest p-value',
    'Gene count',
    paste0('Histogram of stranded vs poly-A RNA-seq\n',
           'differential gene expression RLE normalized\n',
           'RUVg DESeq2 nbinomWaldTest p-values'))

ggsave(
    file.path(
        plot_outdir,
        paste0('stranded_vs_polya_dge_ruvg_k1_',
               'deseq2_nbinom_wald_test_pvals_histogram.png')),
    dpi = 300, plot = p, width = 8, height = 7)


print('Done.')
