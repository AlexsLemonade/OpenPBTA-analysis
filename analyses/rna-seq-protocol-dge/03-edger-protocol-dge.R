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

print('Done.')
