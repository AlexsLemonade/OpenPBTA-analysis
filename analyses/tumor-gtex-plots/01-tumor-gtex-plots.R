# Script to generate tumor vs tumor and tumor vs normal boxplots

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

option_list <- list(
  make_option(c("--expr_mat"), type = "character",
              help = "collapsed TPM expression data: HUGO gene symbol x Sample identifiers (.rds)"),
  make_option(c("--hist_file"), type = "character",
              help = "histologies file (.tsv)"),
  make_option(c("--cohort_list"), type = "character",
              help = "comma separated list of cohorts"),
  make_option(c("--tumor_vs_normal"), type = "logical",
              help = "TRUE or FALSE"),
  make_option(c("--analysis_type"), type = "character",
              help = "cohort_cancer_group_level or cancer_group_level"),
  make_option(c("--plot_width"), type = "character",
              help = "width in pixels"),
  make_option(c("--plot_height"), type = "character",
              help = "height in pixels"),
  make_option(c("--mapping_file"), type = "character",
              help = "filename for writing out file names and other info")
)
# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
expr_mat <- opt$expr_mat
hist_file <- opt$hist_file
cohort_list <- opt$cohort_list
tumor_vs_normal <- opt$tumor_vs_normal
analysis_type <- opt$analysis_type
plot_width <- as.numeric(opt$plot_width)
plot_height <- as.numeric(opt$plot_height)
mapping_file <- opt$mapping_file

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, 'data')
analyses_dir <- file.path(root_dir, 'analyses', 'tumor-gtex-plots')
plots_dir <- file.path(analyses_dir, 'plots')
results_dir <- file.path(analyses_dir, 'results')
dir.create(plots_dir, showWarnings = F, recursive = T)
dir.create(results_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(analyses_dir, 'util', 'pubTheme.R'))
source(file.path(analyses_dir, 'util', 'tumor_plot.R'))
source(file.path(analyses_dir, 'util', 'tumor_vs_normal_plot.R'))

# read input data
expr_mat <- readRDS(expr_mat)
hist_file <- data.table::fread(hist_file)
cohort_list <- trimws(unlist(strsplit(cohort_list,",")))

# filter histologies to cohort_list
hist_file <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% colnames(expr_mat),
         cohort %in% cohort_list) 

# for this module, combine CBTN and PNOC in one cohort = PBTA
hist_file <- hist_file %>%
  mutate(cohort = ifelse(cohort %in% c("CBTN", "PNOC"), "PBTA", cohort))

# match expression matrix to histologies file
expr_mat <- expr_mat %>%
  select(hist_file$Kids_First_Biospecimen_ID)

# just do top 10 genes as a test, we can remove this later
n_genes <- 1
expr_mat <- expr_mat %>%
  head(n = n_genes)

# add gene as a column
expr_mat <- expr_mat %>% 
  rownames_to_column('gene')

# if tumor_vs_normal is TRUE, call tumor_vs_normal_plot else call tumor_plot
# apply function over each gene to generate boxplot and output table
if(tumor_vs_normal){
  print("Tumor vs Normal")
  plyr::d_ply(.data = expr_mat, .variables = "gene", .fun = function(x) tumor_vs_normal_plot(expr_mat_gene = x,
                                                                                             hist_file, analysis_type,
                                                                                             plots_dir, results_dir,
                                                                                             plot_width, plot_height,
                                                                                             mapping_file))
} else {
  print("Tumors only")
  plyr::d_ply(.data = expr_mat, .variables = "gene", .fun = function(x) tumor_plot(expr_mat_gene = x,
                                                                                   hist_file, analysis_type,
                                                                                   plots_dir, results_dir,
                                                                                   plot_width, plot_height,
                                                                                   mapping_file))
}

