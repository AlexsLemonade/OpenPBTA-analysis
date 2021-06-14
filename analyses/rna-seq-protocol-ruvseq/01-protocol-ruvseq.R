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


#------------ Parse parameters -----------------------------
option_list <- list(
    optparse::make_option(
        c("-d", "--dataset"), type = "character",
        help = paste0("Dataset for running differential gene ",
                      "expression analysis: match, dipg, and nbl. "))
)

# parse the parameters
option_parser <- optparse::OptionParser(option_list = option_list)
parsed_opts <- optparse::parse_args(option_parser)
dge_dataset <- parsed_opts$dataset

if (is.null(dge_dataset)) {
    print("Required dataset parameter not found.")
    optparse::print_help(option_parser)
    stop()
}

if (!dge_dataset %in% c('match', 'dipg', 'nbl')) {
    print(paste('Unknown dataset', dge_dataset))
    optparse::print_help(option_parser)
    stop()
}

#------------ Read and pre-process data ----------------------------------------
print('Read count matrices...')

# merged histology df
m_kfbid_sid_htl_df <- readRDS(
    '../../scratch/pbta_kf_gtex_target_tcga_histology_df.rds')

# Kids_First_Biospecimen_ID to read count column ID mapping
kfbid_cid_df <- m_kfbid_sid_htl_df[, c('sample_barcode', 'sample_id')]
colnames(kfbid_cid_df) <- c('Kids_First_Biospecimen_ID', 'col_id')
stopifnot(identical(sum(is.na(kfbid_cid_df)),
                    as.integer(0)))

# load read count matrix and
# select RNA-seq libraries according to @jharenza's comment at
# <https://github.com/PediatricOpenTargets/ticket-tracker/issues/39
#      issuecomment-859751927>
htl_df <- read.delim('../../data/histologies.tsv',
                     stringsAsFactors = FALSE, sep = '\t',
                     header=TRUE)

if (identical(dge_dataset, 'dipg')) {
    cnt_df <- readRDS(
        '../../data/gene-counts-rsem-expected_count-collapsed.rds')
    selected_htl_df <- htl_df %>%
        filter(
            (pathology_diagnosis ==
                 "Brainstem glioma- Diffuse intrinsic pontine glioma") &
                (experimental_strategy == "RNA-Seq") &
                (!sample_id %in% c("7316-1455", "7316-161", "7316-255",
                                   "7316-536", "A16915", "A18777"))
        )
    
    output_dataset_str <- 'dipg_rm_matched_sample_ids'
} else if (identical(dge_dataset, 'match')) {
    cnt_df <- readRDS(
        '../../data/gene-counts-rsem-expected_count-collapsed.rds')
    selected_htl_df <- htl_df %>%
        filter(
            (experimental_strategy == "RNA-Seq") &
                (sample_id %in% c("7316-1455", "7316-161", "7316-255",
                                  "7316-536", "A16915", "A18777"))
        )
    
    output_dataset_str <- 'matched_sample_ids'
} else if (identical(dge_dataset, 'nbl')) {
    cnt_df <- readRDS(
        '../../scratch/pbta_kf_gtex_target_tcga_rsem_expected_cnt_df.rds')
    selected_htl_df <- htl_df %>%
        filter(broad_histology == "Neuroblastoma" &
                   experimental_strategy == "RNA-Seq")
    
    output_dataset_str <- 'nbl'
} else {
    stop(paste0('unknown dataset', dge_dataset))
}
# Find corresponding read count matrix columns of the selected
# RNA-seq libraries
selected_htl_df <- data.frame(
    selected_htl_df[, c('Kids_First_Biospecimen_ID', 'sample_id',
                        'RNA_library')])
selected_htl_df <- merge(selected_htl_df, kfbid_cid_df,
                         by = 'Kids_First_Biospecimen_ID',
                         all.x = TRUE)
stopifnot(identical(sum(is.na(selected_htl_df)),
                    as.integer(0)))
selected_htl_df <- selected_htl_df[order(selected_htl_df$sample_id), ]
rownames(selected_htl_df) <- NULL
# assert
# - there are >= 1 samples in polya or standed
# - RNA_library has only stranded and poly-A values
stopifnot(identical(sort(unique(selected_htl_df$RNA_library)),
                    c("poly-A", "stranded")))
# Subset read count matrix for polya vs stranded DGE analysis
stranded_col_ids <- selected_htl_df[
    selected_htl_df$RNA_library == 'stranded', 'col_id']
polya_col_ids <- selected_htl_df[
    selected_htl_df$RNA_library == 'poly-A', 'col_id']

counts <- cbind(
    cnt_df[, stranded_col_ids],
    cnt_df[, polya_col_ids])

group <- factor(c(rep("stranded", length(stranded_col_ids)),
                  rep("polya", length(polya_col_ids))))

counts_object <- edgeR::DGEList(counts = counts, group = group)
counts_object <- counts_object[
    edgeR::filterByExpr(counts_object), , keep.lib.sizes=FALSE]


#------------ Create output directories ----------------------------------------
table_outdir <- file.path('results', output_dataset_str)
dir.create(table_outdir, showWarnings = FALSE)

plot_outdir <- file.path('plots', output_dataset_str)
dir.create(plot_outdir, showWarnings = FALSE)

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

emp_neg_ctrl_genes <- seg_df$gene[seg_df$gene %in% rownames(round_cnt_mat)]

ruvg_res <- RUVg(seq_expr_set, emp_neg_ctrl_genes, k=1)

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
