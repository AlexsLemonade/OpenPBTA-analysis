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
  make_option(c("--map_file"), type = "character",
              help = "gene symbol-ensembl id mapping file"),
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
  make_option(c("--meta_file"), type = "character",
              help = "filename for writing out file names and other info")
)
# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
expr_mat <- opt$expr_mat
hist_file <- opt$hist_file
map_file <- opt$map_file
cohort_list <- opt$cohort_list
tumor_vs_normal <- opt$tumor_vs_normal
analysis_type <- opt$analysis_type
plot_width <- as.numeric(opt$plot_width)
plot_height <- as.numeric(opt$plot_height)
meta_file <- opt$meta_file

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
source(file.path(analyses_dir, 'util', 'pan_cancer_plot.R'))
source(file.path(analyses_dir, 'util', 'tumor_normal_gtex_plot.R'))

# read input data
expr_mat <- readRDS(expr_mat)
hist_file <- data.table::fread(hist_file)
map_file <- data.table::fread(map_file)
cohort_list <- trimws(unlist(strsplit(cohort_list,",")))

# filter histologies to cohort_list
hist_file <- hist_file %>%
  filter(experimental_strategy == "RNA-Seq",
         cohort == "GTEx" | !is.na(cancer_group), # this filter is needed because TARGET has many normal samples
         Kids_First_Biospecimen_ID %in% colnames(expr_mat),
         cohort %in% cohort_list) 

# match expression matrix to histologies file
expr_mat <- expr_mat %>%
  select(hist_file$Kids_First_Biospecimen_ID)

# just do GPC2 and MYCN as a test, we can remove this later
# expr_mat <- expr_mat[grep('^GPC2$|^MYCN$', rownames(expr_mat)),]

# add gene as a column
expr_mat <- expr_mat %>% 
  rownames_to_column('gene')

# if tumor_vs_normal is TRUE, call tumor_vs_normal_plot else call tumor_plot
# apply function over each gene to generate boxplot and output table
# for-loop performs MUCH faster than plyr::d_ply or other apply functions with large datasets
if(tumor_vs_normal){
  print("Tumor-Normal-GTEx plots")
  for(i in 1:nrow(expr_mat)){
    x <- expr_mat[i,]
    tumor_normal_gtex_plot(expr_mat_gene = x,
                           hist_file, 
                           map_file, 
                           analysis_type,
                           plots_dir, results_dir,
                           plot_width, plot_height,
                           meta_file)
  }
} else {
  print("Pan-cancer plots")
  for(i in 1:nrow(expr_mat)){
    x <- expr_mat[i,]
    pan_cancer_plot(expr_mat_gene = x,
                    hist_file, 
                    map_file,
                    analysis_type,
                    plots_dir, results_dir,
                    plot_width, plot_height,
                    meta_file)
  }
}

