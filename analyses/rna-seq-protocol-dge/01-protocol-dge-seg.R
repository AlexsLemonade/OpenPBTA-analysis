suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))


# functions for upper quartile (UQ) normalization function
UQ <- function(X){
    uq <- function(y){
        quantile(y, 0.75)
    }

    X <- X + 0.1
    upperQ <- apply(X, 2, uq)
    f <- upperQ / mean(upperQ) # calculate normalization factor
    res <- scale(X, center=FALSE, scale=f) 
    return(res)
}
# per gene normalization by Median: pgQ2
# X: a matrix of data with the multiplication of factor (f) as:
# f=100 as a default
uq.pgQ2 <- function(X, f = 100){
    uq.res<-UQ(X) #  perform UQ normalization per sample

    ## perform per gene nomralization
    m<-apply(uq.res, 1, median)
    
    idx <- which(m == 0) # avoid 0 median 
    m[idx] = 0.1
    
    si <- m/f # calculate normalization factor
    X1 <- scale(t(uq.res), center=FALSE, scale=si)
    res <- t(X1)
    return(res)  
}
#------------ Parse parameters -----------------------------
option_list <- list(
  make_option(c("-n", "--normalization"), type = "character",
              help = paste0("RSEM expected count normalization method: ",
                            "tmm or uqpgq2. "))
)

# parse the parameters
option_parser <- OptionParser(option_list = option_list)
parsed_opts <- parse_args(option_parser)
norm_method <- parsed_opts$normalization

if (is.null(norm_method)) {
    print("Required normalization parameter not found.")
    print_help(option_parser)
    stop()
}

if (!norm_method %in% c('tmm', 'uqpgq2')) {
    print(paste('Unknown normalization method', norm_method))
    print_help(option_parser)
    stop()
}
#------------ Create output directories --------------------
table_outdir <- file.path('results', paste0(norm_method, '_normalized'))
dir.create(table_outdir, showWarnings = FALSE)

plot_outdir <- file.path('plots', paste0(norm_method, '_normalized'))
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
# Make design matrix 
design <- model.matrix(~group)

# sample_prot is used to select stably expressed genes
sample_prot <- c('s1_stranded', 's2_stranded', 's3_stranded',
                 's4_stranded', 's5_stranded', 's6_stranded',
                 's1_polya', 's2_polya', 's3_polya',
                 's4_polya', 's5_polya', 's6_polya')

### build normalized DGEList
counts_object <- DGEList(counts = counts, group = group)
counts_object <- counts_object[
    filterByExpr(counts_object), , keep.lib.sizes=FALSE]
# read housekeeping genes
hkgenes <- read.csv(file.path('input', 'Housekeeping_GenesHuman.csv'),
                    sep = ";", stringsAsFactors = FALSE)
hkgenes <- unique(hkgenes$Gene.name)
hkgenes <- intersect(rownames(counts_object), hkgenes)

if (identical(norm_method, 'tmm')) {
    counts_object <- calcNormFactors(counts_object, method="TMM")
    norm_cnt_mat <- cpm(counts_object)

    counts_object <- counts_object[hkgenes, ]
    norm_cnt_mat <- norm_cnt_mat[hkgenes, ]
    round_norm_cnt_mat <- round(norm_cnt_mat)

    plot_nm_str <- 'TMM'
} else if (identical(norm_method, 'uqpgq2')) {
    counts_object <- DGEList(counts = uq.pgQ2(X = counts_object$counts),
                             group = group)
    counts_object <- counts_object[hkgenes, ]

    norm_cnt_mat <- counts_object$counts
    round_norm_cnt_mat <- round(norm_cnt_mat)
    counts_object <- calcNormFactors(counts_object, method="none")

    plot_nm_str <- 'UQ-pgQ2'
}
### Estimate dispersion
counts_object <- estimateDisp(
    counts_object, design=design, prior.df=NULL, trend.method='locfit',
    tagwise=TRUE)
png(file.path(plot_outdir, 'estimated_dispersions.png'))
plotBCV(counts_object)
n_null_dev <- dev.off()
#------------ Run LRT --------------------
print(paste0('Run differential gene expression LRT on',
             ' poly-A vs stranded RNA-seq...'))
glm_fit <- glmFit(counts_object, design)
lrt <- glmLRT(glm_fit)

# Save LRT p-value table and histogram
lrt_top_tags <- topTags(lrt, adjust.method = "BH", n = Inf)
lrt_out_df <- lrt_top_tags$table
stopifnot(identical(colnames(lrt_out_df),
                    c("logFC", "logCPM", "LR", "PValue", "FDR")))
colnames(lrt_out_df) <- c("stranded_over_polya_logFC",
                          "average_logCPM", "LR", "PValue", "BH FDR")
write.csv(lrt_out_df, file=file.path(table_outdir,
                                     'stranded_vs_polya_dge_lrt_res.csv'))

# plot and save p-value histogram
p <- ggplot(lrt_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab(paste0('stranded vs poly-A RNA-seq DGE ', plot_nm_str,
                ' LRT p-value')) +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression ', plot_nm_str,
                   ' normalized\n',
                   'LRT p-values\n',
                   sum(lrt_top_tags$table$PValue < 0.05),
                   ' genes have p-value < 0.05\n',
                   sum(lrt_top_tags$table$PValue >= 0.05),
                   ' genes have p-value >= 0.05\n',
                   sum(lrt_top_tags$table$FDR < 0.05),
                   ' genes have BH FDR < 0.05\n',
                   sum(lrt_top_tags$table$FDR >= 0.05),
                   ' genes have BH FDR >= 0.05')) +
    theme(text = element_text(size=15))

ggsave(file.path(plot_outdir, 'stranded_vs_polya_dge_lrt_pvals_histogram.png'),
       dpi = 300, plot = p, width = 8, height = 7)
#------------ Run exactTest --------------------
print(paste0('Run differential gene expression exactTest on',
             ' poly-A vs stranded RNA-seq...'))
et <- exactTest(counts_object)
# get p-value table
et_top_tags <- topTags(et, adjust.method = "BH", n = Inf)
et_out_df <- et_top_tags$table
stopifnot(identical(colnames(et_out_df),
                    c("logFC", "logCPM", "PValue", "FDR")))
colnames(et_out_df) <- c("stranded_over_polya_logFC",
                         "average_logCPM", "PValue", "BH FDR")
write.csv(et_out_df,
          file=file.path(table_outdir,
                         'stranded_vs_polya_dge_exact_test_res.csv'))

# plot and save p-value histogram
p <- ggplot(et_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab(paste0('stranded vs poly-A RNA-seq DGE ', plot_nm_str,
                ' exactTest p-value')) +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression ', plot_nm_str,
                   ' normalized\n',
                   'exactTest p-values\n',
                   sum(et_top_tags$table$PValue < 0.05),
                   ' genes have p-value < 0.05\n',
                   sum(et_top_tags$table$PValue >= 0.05),
                   ' genes have p-value >= 0.05\n',
                   sum(et_top_tags$table$FDR < 0.05),
                   ' genes have BH FDR < 0.05\n',
                   sum(et_top_tags$table$FDR >= 0.05),
                   ' genes have BH FDR >= 0.05')) +
    theme(text = element_text(size=15))

ggsave(file.path(plot_outdir,
                 'stranded_vs_polya_dge_exact_test_pvals_histogram.png'),
       dpi = 300, plot = p, width = 8, height = 7)
#------------ Run QL F test ------------------------------------
print(paste0('Run differential gene expression QL F-test on',
             ' poly-A vs stranded RNA-seq...'))
ql_fit <- glmQLFit(counts_object, design)
png(file.path(plot_outdir, 'glm_ql_fit_deviances.png'))
plotQLDisp(ql_fit)
n_null_dev <- dev.off()
qlft <- glmQLFTest(ql_fit)
# get p-value table
qlft_top_tags <- topTags(qlft, adjust.method = "BH", n = Inf)
qlft_out_df <- qlft_top_tags$table
stopifnot(identical(colnames(qlft_out_df),
                    c("logFC", "logCPM", "F", "PValue", "FDR")))
colnames(qlft_out_df) <- c("stranded_over_polya_logFC",
                           "average_logCPM", "F", "PValue", "BH FDR")
write.csv(qlft_out_df,
          file=file.path(table_outdir,
                         'stranded_vs_polya_dge_ql_ftest_res.csv'))
# plot and save p-value histogram
p <- ggplot(qlft_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab(paste0('stranded vs poly-A RNA-seq DGE ', plot_nm_str,
                ' QL F-test p-value')) +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression ', plot_nm_str,
                   ' normalized\n',
                   'QL F-test p-values\n',
                   sum(qlft_top_tags$table$PValue < 0.05),
                   ' genes have p-value < 0.05\n',
                   sum(qlft_top_tags$table$PValue >= 0.05),
                   ' genes have p-value >= 0.05\n',
                   sum(qlft_top_tags$table$FDR < 0.05),
                   ' genes have BH FDR < 0.05\n',
                   sum(qlft_top_tags$table$FDR >= 0.05),
                   ' genes have BH FDR >= 0.05')) +
    theme(text = element_text(size=15))

ggsave(file.path(plot_outdir,
                 'stranded_vs_polya_dge_ql_ftest_pvals_histogram.png'),
       dpi = 300, plot = p, width = 8, height = 7)
#------------ Run DESeq2 nbinomWaldTest ------------------------
print(paste0('Run differential gene expression DESeq2 nbinomWaldTest on',
             ' poly-A vs stranded RNA-seq...'))
suppressMessages(
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round_norm_cnt_mat,
        colData = data.frame(group),
        design = ~ group)
)
DESeq2::sizeFactors(dds) <- rep(1, dim(round_norm_cnt_mat)[2])
suppressMessages(dds <- DESeq2::estimateDispersions(dds, fitType='local'))
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
    xlab(paste0('stranded vs poly-A RNA-seq DGE ', plot_nm_str,
                ' nbinomWaldTest p-value')) +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression ', plot_nm_str,
                   ' normalized\n',
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
#------------ Select stably expressed genes --------------------
print('Select stably expressed genes in stranded and poly-A libraries...')
# merge exact test result, LRT result, QL F-test result, and read counts into
# a single data_frame

# exact test data frame
et_m_df <- et_out_df[, c('stranded_over_polya_logFC',
                         'average_logCPM', 'PValue')]
colnames(et_m_df) <- c('stranded_over_polya_logFC',
                       'average_logCPM', 'exact_test_pval')
et_m_df <- rownames_to_column(et_m_df, 'gene')
# LRT data frame
lrt_m_df <- lrt_out_df[, c('PValue'), drop=FALSE]
colnames(lrt_m_df) <- c('lrt_pval')
lrt_m_df <- rownames_to_column(lrt_m_df, 'gene')
# QL F-test data frame
qlft_m_df <- qlft_out_df[, c('PValue'), drop=FALSE]
colnames(qlft_m_df) <- c('qlft_pval')
qlft_m_df <- rownames_to_column(qlft_m_df, 'gene')
# DESeq2 nbinomWaldTest
deseq2_nbwt_m_df <- deseq2_out_df[, c('PValue'), drop=FALSE]
colnames(deseq2_nbwt_m_df) <- c('deseq2_nbwaldtest_pval')
deseq2_nbwt_m_df <- rownames_to_column(deseq2_nbwt_m_df, 'gene')
# some DESeq2 p-values are NAs due to all-0 normalized read counts
deseq2_nbwt_m_df[is.na(deseq2_nbwt_m_df)] <- 1

# normalized read count data frame
norm_cnt_df <- data.frame(norm_cnt_mat)
cnt_m_df <- norm_cnt_df[lrt_m_df$gene, ]
# assert samples in cnt_m_df are in the same order as sample_prot
stopifnot(identical(colnames(cnt_m_df),
                    c("BS_HE0WJRW6", "BS_SHJA4MR0", "BS_FN07P04C",
                      "BS_8QB4S4VA", "BS_KABQQA0T", "BS_D7XRFE0R",
                      "BS_HWGWYCY7", "BS_X0XXN9BK", "BS_W4H1D4Y6",
                      "BS_QKT3TJVK", "BS_7WM3MNZ0", "BS_68KX6A42")))
# sample_prot defined at top when reading data
colnames(cnt_m_df) <- paste(sample_prot, colnames(norm_cnt_df),
                            sep='_')
# coefficient of variation data frame
coeff_var <- function(x, na.rm = FALSE)  {
    return(sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
}
coeff_var_df <- data.frame(normalized_count_cv = apply(cnt_m_df, 1, coeff_var))
coeff_var_df <- rownames_to_column(coeff_var_df, 'gene')

cnt_m_df <- rownames_to_column(cnt_m_df, 'gene')
# merge
seg_df <- Reduce(function(x, y) merge(x, y, by = 'gene'),
                 list(et_m_df, lrt_m_df, qlft_m_df, deseq2_nbwt_m_df,
                      coeff_var_df, cnt_m_df))

# plot and save CV histogram
p <- ggplot(seg_df, aes(x=normalized_count_cv)) +
    geom_histogram(binwidth = 0.1, center = 0.05) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02)),
                       breaks = seq(0, 3, 0.5)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab('Normalized count coefficient of variation\nof each gene') +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of normalized count coefficient\n',
                   'of variation of each gene')) +
    theme(text = element_text(size=15))
ggsave(file.path(plot_outdir, 'normalized_count_gene_cv_histogram.png'),
       dpi = 300, plot = p, width = 8, height = 5)

# plot and save logFC histogram
p <- ggplot(seg_df, aes(x=stranded_over_polya_logFC)) +
    geom_histogram(binwidth = 0.2, center = 0.1) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02)),
                       limits = c(-8, 8)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab('edgeR DGE logFC of each gene') +
    ylab('Gene count') +
    ggtitle('Histogram of edgeR DGE logFC ') +
    theme(text = element_text(size=15))
suppressWarnings(ggsave(file.path(plot_outdir,
                                  'edger_logfc_histogram.png'),
                        dpi = 300, plot = p, width = 8, height = 5))

seg_df <- seg_df[
    seg_df$normalized_count_cv < quantile(seg_df$normalized_count_cv, 0.25), ]
# select all pval > 0.05
seg_df <- seg_df[seg_df$lrt_pval > 0.05 &
                     seg_df$exact_test_pval > 0.05 &
                     seg_df$qlft_pval > 0.05 &
                     seg_df$deseq2_nbwaldtest_pval > 0.05, ]
seg_df <- seg_df[abs(seg_df$stranded_over_polya_logFC) <
                     median(abs(seg_df$stranded_over_polya_logFC)), ]
# select mean expression between 25th and 75th percentile
alc_idc <- seg_df$average_logCPM > quantile(seg_df$average_logCPM, 0.25)
alc_idc <- alc_idc & (seg_df$average_logCPM < 
                          quantile(seg_df$average_logCPM, 0.75))
seg_df <- seg_df[alc_idc, ]
# order by CV
seg_df <- seg_df[order(seg_df$normalized_count_cv), ]
rownames(seg_df) <- NULL

write.csv(seg_df,
          file = file.path(table_outdir,
                           'stranded_vs_polya_stably_exp_genes.csv'))

get_ge_boxplot <- function(df_row) {
    stopifnot(identical(nrow(df_row), as.integer(1)))
    gid <- df_row[, 'gene']
    plot_df <- df_row[, c("s1_stranded_BS_HE0WJRW6", 
                         "s2_stranded_BS_SHJA4MR0", "s3_stranded_BS_FN07P04C", 
                         "s4_stranded_BS_8QB4S4VA", 
                         "s5_stranded_BS_KABQQA0T", "s6_stranded_BS_D7XRFE0R",
                         "s1_polya_BS_HWGWYCY7", 
                         "s2_polya_BS_X0XXN9BK", "s3_polya_BS_W4H1D4Y6", 
                         "s4_polya_BS_QKT3TJVK", 
                         "s5_polya_BS_7WM3MNZ0", "s6_polya_BS_68KX6A42")]
    # wide to long
    plot_df <- gather(plot_df)
    sid_match_res <- str_match(plot_df$key, '^(s[0-9])_([a-z]+)_(.+)$')
    plot_df$sample_id <- sid_match_res[, 2]
    plot_df$protocol <- sid_match_res[, 3]
    plot_df$bio_id <- sid_match_res[, 4]

    p <- plot_df %>%
        ggplot(aes(protocol, value, label = bio_id)) +
        geom_boxplot() +
        geom_point() +
        geom_line(aes(group=sample_id)) +
        geom_text(check_overlap = FALSE, angle = 0, hjust = 0,
                  nudge_x = 0.03, size=3) +
        theme(legend.position = "none") +
        xlab('RNA-seq library protocol') +
        ylab('Normalized read count') +
        ggtitle(paste0(gid, '\nLines connect biospecimens\n',
                       'with the same sample IDs')) +
        theme_classic() +
        theme(text = element_text(size=15))
    return(p)
}

seg_plot_outdir <- file.path(plot_outdir,
                             'stably_exp_gene_protocol_diff_boxplot')
dir.create(seg_plot_outdir, showWarnings = FALSE)

res <- sapply(1:30, function(i) {
    gid <- seg_df[i, 'gene']
    p <- get_ge_boxplot(seg_df[i, ])
    ggsave(file.path(seg_plot_outdir,
                     paste0(gid, 
                            '_normalized_count_protocol_diff_boxplot.png')),
           dpi = 300, plot = p, width = 5, height = 5)
    return(TRUE)
})
print('Done.')
