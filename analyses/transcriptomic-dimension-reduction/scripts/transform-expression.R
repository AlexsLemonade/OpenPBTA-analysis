# Taroni for CCDL 2019
# Takes an expression matrix and writes a log2(x + 1) transformed matrix to file
#
# USAGE:
#   Rscript --vanilla analyses/transcriptomic-dimension-reduction/scripts/transform-expression.R \
#     --expression data/pbta-gene-expression-kallisto.stranded.rds \
#     --output scratch/log2-kallisto.stranded.rds

`%>%` <- dplyr::`%>%`

option_list <- list(
  optparse::make_option(
    c("-i", "--expression"),
    type = "character",
    default = NULL,
    help = "full path to expression data RDS file",
  ),
  optparse::make_option(
    c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "output file"
  )
)

expression_file <- opt$expression
output_file <- opt$output_file

expression_data <- readr::read_rds(expression_file) %>%
    as.data.frame()

# the first column will be either transcript or column ids
feature_identifier <- colnames(expression_data)[1]

# drop any columns that contain other identifers
expression_matrix <- expression_data %>%
  dplyr::select(!!rlang::sym(feature_identifier), dplyr::starts_with("BS_")) %>%
  tibble::column_to_rownames(var = feature_identifier) %>%
  as.matrix()

# the transform step itself
log2_expression_matrix <- log2(expression_matrix + 1)

# include all feature identifiers in what gets written to file
log2_expression_df <- as.data.frame(
  cbind(expression_data[, grepl("transcript_id|gene_id",
                                colnames(expression_data))],
        log2_expression_matrix)
)

# write to file
readr::write_rds(log2_expression_df, output_file)
