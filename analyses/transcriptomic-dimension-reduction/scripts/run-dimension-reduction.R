# Unsupervised Analysis of Transcriptomic Differences - Run Dimension Reduction
# Chante Bethell for CCDL 2019
#
# This script runs dimensionality reduction techniques:
# Principal Component Analysis (PCA), Uniform Manifold Approximation
# and Projection (UMAP), and optionally t-Distributed Stochastic Neighbor
# Embedding (t-SNE). It will perform t-SNE by default.
#
# It addresses issue #9 (https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/9)
# in the Open-PBTA analysis repository.
#
# Usage
#
# This script is intended to be run via the command line from the top directory
# of the repository as follows (for all RSEM samples, with t-SNE skipped):
#
# Rscript --vanilla analyses/transcriptomic-dimension-reduction/scripts/run-dimension-reduction.R \
#   --expression data/pbta-gene-expression-rsem.fpkm.rds \
#   --metadata data/pbta-histologies.tsv \
#   --filename_lead rsem_all \
#   --output_directory analyses/transcriptomic-dimension-reduction/results \
#   --seed 2019 \
#   --perplexity 10 \
#   --neighbors 15 \
#   --low_count_threshold 100 \
#   --skip_tsne

#### Install packages ----------------------------------------------------------

# This is needed for running the t-SNE analysis
if (!("Rtsne" %in% installed.packages())) {
  install.packages("Rtsne")
}

# This is needed for running the umap analysis
if (!("umap" %in% installed.packages())) {
  install.packages("umap")
}

# This is needed for taking arguments from the command line
if (!("optparse" %in% installed.packages())) {
  install.packages("optparse")
}

#### Functions -----------------------------------------------------------------

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# source the functions for dimension reduction
source(file.path(root_dir, "analyses", "transcriptomic-dimension-reduction",
                 "util", "dimension-reduction-functions.R"))

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--expression"),
    type = "character",
    default = NULL,
    help = "full path to expression data RDS file",
  ),
  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to metadata TSV file"
  ),
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory"
  ),
  optparse::make_option(
    c("-f", "--filename_lead"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2019,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-p", "--perplexity"),
    type = "integer",
    default = 10,
    help = "t-SNE perplexity parameter integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-n", "--neighbors"),
    type = "integer",
    default = 15,
    help = "UMAP neighbors parameter integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-l", "--low_count_threshold"),
    type = "integer",
    default = 100,
    help = "Number of total counts required to retain a gene",
    metavar = "integer"
  ),
  optparse::make_option(
    "--skip_tsne",
    action = "store_true",
    default = FALSE,
    help = "If used, t-SNE will be skipped"
  ),
  optparse::make_option(
    "--log2_transform",
    action = "store_true",
    default = FALSE,
    help = "If used, log2(x + 1) transformation will be performed prior to dimension reduction"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# RDS file that contains the expression data
expression_file <- opt$expression

# TSV file that contains metadata
metadata_file <- opt$metadata

# Create specified output directory if it does not yet exist
output_directory <- opt$output_directory
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Character string that will be used to name files
filename_lead <- opt$filename_lead

# Define seed as the parameter passed via the command line
seed <- opt$seed

# Define t-SNE perplexity parameter as the parameter passed via command line
perplexity_parameter <- opt$perplexity

# Define UMAP neighbors parameter as the parameter passed via command line
neighbors_parameter <- opt$neighbors

# Threshold used to filter out genes with low counts
gene_count_threshold <- opt$low_count_threshold

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata_df <- data.frame(readr::read_tsv(metadata_file))

# Read in expression data
expression_data <- readr::read_rds(expression_file) %>%
  as.data.frame()

# the first column will be either transcript or column ids
feature_identifier <- colnames(expression_data)[1]

# drop any columns that contain other identifers
expression_data <- expression_data %>%
  dplyr::select(!!rlang::sym(feature_identifier), dplyr::starts_with("BS_")) %>%
  dplyr::filter(complete.cases(.)) %>%
  tibble::column_to_rownames(var = feature_identifier)

# Filter out low count genes
genes_to_keep <- rowSums(expression_data) >= gene_count_threshold
expression_data <- expression_data[genes_to_keep, ]

if (opt$log2_transform) {
  expression_data <- log2(expression_data + 1)
}

# Transpose the data
transposed_exp_data <- t(expression_data)

#### Dimension reduction -------------------------------------------------------

# PCA
pca_file <- paste0(filename_lead, "_pca")
dimension_reduction_wrapper(
  transposed_expression_matrix = transposed_exp_data,
  method = "PCA",
  seed = seed,
  metadata_df = metadata_df,
  filename_lead = pca_file,
  output_directory = output_directory
)

# UMAP
umap_file <- paste0(filename_lead, "_umap")
dimension_reduction_wrapper(
  transposed_expression_matrix = transposed_exp_data,
  method = "UMAP",
  seed = seed,
  metadata_df = metadata_df,
  filename_lead = umap_file,
  output_directory = output_directory,
  neighbors_parameter = neighbors_parameter
)

# t-SNE, if not skipped
if (!opt$skip_tsne) {
  tsne_file <- paste0(filename_lead, "_tsne")
  dimension_reduction_wrapper(
    transposed_expression_matrix = transposed_exp_data,
    method = "t-SNE",
    seed = seed,
    metadata_df = metadata_df,
    filename_lead = tsne_file,
    output_directory = output_directory,
    perplexity_parameter = perplexity_parameter
  )
}
