suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))

#------------ Read and pre-process data --------------------
print('Read count matrices...')
polya <- readRDS(
    "results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds")
stranded <- readRDS(
    "results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds")

common_gene_id <- intersect(rownames(stranded), rownames(polya))
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
    stranded[common_gene_id, c("BS_HE0WJRW6", "BS_SHJA4MR0", "BS_FN07P04C",
                               "BS_8QB4S4VA", "BS_KABQQA0T", "BS_D7XRFE0R")],
    polya[common_gene_id, c("BS_HWGWYCY7", "BS_X0XXN9BK", "BS_W4H1D4Y6",
                            "BS_QKT3TJVK", "BS_7WM3MNZ0", "BS_68KX6A42")]
)
group <- factor(c(rep("stranded", 6), rep("polya", 6)))

# sample_prot is used to select stably expressed genes
sample_prot <- c('s1_stranded', 's2_stranded', 's3_stranded',
                 's4_stranded', 's5_stranded', 's6_stranded',
                 's1_polya', 's2_polya', 's3_polya',
                 's4_polya', 's5_polya', 's6_polya')

### build uq.pgQ2 normalized DGEList
counts_object = DGEList(counts = counts, group = group)
counts_object <- counts_object[
    filterByExpr(counts_object), , keep.lib.sizes=FALSE]
# functions for upper quartile (UQ) normalization function
# copied from
# <https://github.com/d3b-center/OMPARE/blob/
#  2c5cb9f56dc95212f113a6d13459e86b9d70c48f/code/
#  patient_level_analyses/utils/rnaseq_edger_normalizations.R#L3-L33>
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
counts_object$counts <- uq.pgQ2(X = counts_object$counts)

### Estimate common dispersion from a set of housekeeping genes.
hkgenes <- read.csv(file.path('input', 'Housekeeping_GenesHuman.csv'),
                    sep = ";", stringsAsFactors = FALSE)
hkgenes <- unique(hkgenes$Gene.name)
hkgenes <- intersect(rownames(counts_object), hkgenes)
hkg_counts_object <- counts_object[hkgenes, ]
# the expression prints "Design matrix not provided. Switch to the classic 
# mode.". Suppress this message when running.
out_msg <- capture.output(
    hkg_counts_object <- estimateDisp(hkg_counts_object, trend.method="none",
                                      tagwise=FALSE))
counts_object$common.dispersion <- hkg_counts_object$common.dispersion
### Make design matrix 
design <- model.matrix(~group)
#------------ Run LRT --------------------
print(paste0('Run differential gene expression LRT on',
             ' poly-A vs stranded RNA-seq...'))
fit <- glmFit(counts_object, design)
lrt <- glmLRT(fit)

# Save LRT p-value table and histogram
lrt_top_tags <- topTags(lrt, adjust.method = "BH", n = Inf)
lrt_out_df <- lrt_top_tags$table
stopifnot(identical(colnames(lrt_out_df),
                    c("logFC", "logCPM", "LR", "PValue", "FDR")))
colnames(lrt_out_df) <- c("stranded_over_polya_logFC",
                          "average_logCPM", "LR", "PValue", "BH FDR")
write.csv(lrt_out_df, file='results/stranded_vs_polya_dge_lrt_res.csv')

# plot and save p-value histogram
p <- ggplot(lrt_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab('stranded vs poly-A RNA-seq DGE UQ-pgQ2 LRT p-value') +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression UQ-pgQ2 normalized\n',
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

ggsave('plots/stranded_vs_polya_dge_lrt_pvals_histogram.png', dpi = 300,
       plot = p, width = 8, height = 7)


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
write.csv(et_out_df, file='results/stranded_vs_polya_dge_exact_test_res.csv')

# plot and save p-value histogram
p <- ggplot(et_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab('stranded vs poly-A RNA-seq DGE UQ-pgQ2 exactTest p-value') +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression UQ-pgQ2 normalized\n',
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

ggsave('plots/stranded_vs_polya_dge_exact_test_pvals_histogram.png', dpi = 300,
       plot = p, width = 8, height = 7)
#------------ Run QL F test ------------------------------------
print(paste0('Run differential gene expression QL F-test on',
             ' poly-A vs stranded RNA-seq...'))
ql_fit <- glmQLFit(counts_object, design)
qlft <- glmQLFTest(ql_fit)
# get p-value table
qlft_top_tags <- topTags(qlft, adjust.method = "BH", n = Inf)
qlft_out_df <- qlft_top_tags$table
stopifnot(identical(colnames(qlft_out_df),
                    c("logFC", "logCPM", "F", "PValue", "FDR")))
colnames(qlft_out_df) <- c("stranded_over_polya_logFC",
                           "average_logCPM", "F", "PValue", "BH FDR")
write.csv(qlft_out_df,
          file='results/stranded_vs_polya_dge_ql_ftest_res.csv')
# plot and save p-value histogram
p <- ggplot(qlft_top_tags$table, aes(x=PValue)) +
    geom_histogram(binwidth = 0.05, center = 0.025) +
    theme_classic() +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    xlab('stranded vs poly-A RNA-seq DGE UQ-pgQ2 QL F-test p-value') +
    ylab('Gene count') +
    ggtitle(paste0('Histogram of stranded vs poly-A RNA-seq\n',
                   'differential gene expression UQ-pgQ2 normalized\n',
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

ggsave('plots/stranded_vs_polya_dge_ql_ftest_pvals_histogram.png', dpi = 300,
       plot = p, width = 8, height = 7)

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

# normalized read count data frame
norm_cnt_df <- data.frame(counts_object$counts)
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
                 list(et_m_df, lrt_m_df, qlft_m_df, coeff_var_df, cnt_m_df))

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
ggsave('plots/normalized_count_gene_cv_histogram.png', dpi = 300,
       plot = p, width = 8, height = 5)

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
suppressWarnings(ggsave('plots/edger_logfc_histogram.png', dpi = 300,
                        plot = p, width = 8, height = 5))

# select CV < 0.5
seg_df <- seg_df[seg_df$normalized_count_cv < 0.5, ]
# select pval > 0.05
seg_df <- seg_df[seg_df$lrt_pval > 0.05 & seg_df$exact_test_pval > 0.05 &
                     seg_df$qlft_pval > 0.05, ]
# select abs(logFC) < 0.5
seg_df <- seg_df[abs(seg_df$stranded_over_polya_logFC) < 0.5, ]
# select average_logCPM > 0.75
seg_df <- seg_df[seg_df$average_logCPM > 0.8, ]
# order by CV
seg_df <- seg_df[order(seg_df$normalized_count_cv), ]
rownames(seg_df) <- NULL

write.csv(seg_df, file = 'results/stranded_vs_polya_stably_exp_genes.csv')

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

seg_plot_outdir <- file.path('plots', 'stably_exp_gene_protocol_diff_boxplot')
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
