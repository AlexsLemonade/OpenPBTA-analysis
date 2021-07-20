# This script adds gene and cancer_group annotations to a long-format table
#
# This script parses arguments and calls the annotate_long_format_table function
# in the annotator/annotator-api.R file to add required annotation columns
#
# Example usages:



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
    c("-l", "--long-format-table-tsv"), type = "character",
    help = paste0("Path to the long-format table TSV file to be annotated"))
)
option_parser <- optparse::OptionParser(option_list = option_list)
parsed_opts <- optparse::parse_args(option_parser)
print(parsed_opts)
