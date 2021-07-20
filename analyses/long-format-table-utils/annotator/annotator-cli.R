# This script adds gene and cancer_group annotations to an input long-format
# table TSV file and outputs an annotated long-format table TSV file
#
# This script parses arguments and calls the annotate_long_format_table function
# in the annotator/annotator-api.R file to add required annotation columns
#
# EXAMPLE USAGES:
#
# - Print help message:
#
# Rscript --vanilla analyses/long-format-table-utils/annotator/annotator-cli.R \
#   -h
#
# - Add RMTL, EFO, and MONDO columns
#   - The `-r` option replaces NAs with empty strings for **ALL** columns of the
#     input table
#   - The `-v` option prints extra messages on progress
#
# Rscript --vanilla analyses/long-format-table-utils/annotator/annotator-cli.R \
#   -r -v -c RMTL,EFO,MONDO \
#   -i long_n_tpm_mean_sd_quantile_group_gene_wise_zscore.tsv \
#   -o long_n_tpm_mean_sd_quantile_group_gene_wise_zscore_annotated.tsv



# source annotator-api.R to get the annotate_long_format_table function --------
# Detect the ".git" folder -- this will in the project root directory. Use
# this as the root directory to ensure proper execution, no matter where it is
# called from.
#
# This only works if the working directory is OpenPedCan-analysis or a
# subdirectory of OpenPedCan-analysis
#
# root_dir is the absolute path of OpenPedCan-analysis
#
# Adapted from the oncoprint-landscape module.
tryCatch(
  {
    root_dir <- rprojroot::find_root(rprojroot::has_file(".git/index"))
  },
  error = function(err_cond) {
    # adapted from http://adv-r.had.co.nz/Exceptions-Debugging.html
    err_cond$message <- paste0(
      err_cond$message,
      "\nTry re-running this script with working directory as ",
      "OpenPedCan-analysis or a subdirectory of OpenPedCan-analysis.\n")
    stop(err_cond)
  }
)
# Get the annotate_long_format_table function from annotator-api.R
source(file.path(
  root_dir, "analyses", "long-format-table-utils", "annotator",
  "annotator-api.R"))



# Parse arguments --------------------------------------------------------------
option_list <- list(
  optparse::make_option(
    c("-r", "--replace-na-with-empty-string"), action = "store_true",
    default = FALSE,
    help = paste0(
      "Replace NAs with empty strings for **ALL** columns of the input table")),
  optparse::make_option(
    c("-c", "--columns-to-add"), type = "character",
    default = paste0(
      "RMTL,Gene_type,OncoKB_cancer_gene,OncoKB_oncogene_TSG,",
      "Gene_full_name,Protein_RefSeq_ID,EFO,MONDO"),
    help = paste0(
      "A comma-separated list of unique annotation columns to be added, ",
      "e.g. \"EFO,MONDO\" and \"RMTL,Gene_type,OncoKB_cancer_gene\". ",
      "Available columns are: RMTL, Gene_type, OncoKB_cancer_gene, ",
      "OncoKB_oncogene_TSG, Gene_full_name, Protein_RefSeq_ID, EFO, MONDO. ",
      "[Default value is \"%default\", which is to add all available ",
      "annotation columns]")),
  optparse::make_option(
    c("-i", "--input-long-format-table-tsv"), type = "character",
    help = "Path to the input long-format table TSV file to be annotated"),
  optparse::make_option(
    c("-o", "--output-long-format-table-tsv"), type = "character",
    help = "Path to output the annotated long-format table TSV file"),
  optparse::make_option(
    c("-v", "--verbose"), action = "store_true",
    default = FALSE,
    help = "Print extra messages on progress")
)
option_parser <- optparse::OptionParser(option_list = option_list)
parsed_opts <- optparse::parse_args(option_parser)



# Read, annotate, and output ---------------------------------------------------
columns_to_add <- stringr::str_split(parsed_opts$`columns-to-add`, ",")[[1]]

if (parsed_opts$verbose) {
  message(paste0("Read ", parsed_opts$`input-long-format-table-tsv`, "..."))
}
input_df <- readr::read_tsv(
  parsed_opts$`input-long-format-table-tsv`, na = character(),
  col_types = readr::cols(.default = readr::col_guess()))

if (parsed_opts$verbose) {
  message(paste0("Annotate ", parsed_opts$`input-long-format-table-tsv`, "..."))
}
ann_df <- annotate_long_format_table(
  long_format_table = input_df, columns_to_add = columns_to_add,
  replace_na_with_empty_string = parsed_opts$`replace-na-with-empty-string`)

if (parsed_opts$verbose) {
  message(paste0("Output ", parsed_opts$`output-long-format-table-tsv`, "..."))
}
readr::write_tsv(ann_df, parsed_opts$`output-long-format-table-tsv`)

if (parsed_opts$verbose) {
  message("Done.")
}
