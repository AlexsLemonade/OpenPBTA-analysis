# Add gene and cancer_group annotations to long-format tables
#
# Args:
# - long_format_table: A data.frame or tibble that contains the following
#   columns as character:
#   - Gene_symbol: HUGO symbols, e.g. PHLPP1, TM6SF1, and DNAH5
#   - Gene_Ensembl_ID: Ensembl ENSG IDs without `.#` versions, e.g.
#     ENSG00000039139, ENSG00000111261, and ENSG00000169710
#   - Disease: The `cancer_group` in the `histologies.tsv`, e.g.
#     Adamantinomatous Craniopharyngioma, Atypical Teratoid Rhabdoid Tumor, and
#     Low-grade glioma/astrocytoma
# - is_gene_level_table: TRUE or FALSE on whether the table is gene-level. If is
#   gene-level, Gene_type column will be added, otherwise will not. Default
#   value is FALSE.
# - add_Protein_RefSeq_ID: TRUE or FALSE on whether to add Protein_RefSeq_ID
#   column. Default value is FALSE.
#
# Returns a data.frame or tibble, based on input table type, with additonal
# annotation columns
annotate_long_format_table <- function(long_format_table,
                                       is_gene_level_table = FALSE,
                                       add_Protein_RefSeq_ID = FALSE) {
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
}
