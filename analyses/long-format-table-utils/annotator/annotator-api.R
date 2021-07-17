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
# - is_gene_level_table: TRUE or FALSE on whether the input table is gene-level.
#   If the input table is gene-level, Gene_type column will be added, otherwise
#   will not. Default value is FALSE.
# - add_Protein_RefSeq_ID: TRUE or FALSE on whether to add Protein_RefSeq_ID
#   column. Default value is FALSE.
# - replace_na_with_empty_string: TRUE or FALSE on whether to replace NAs with
#   empty strings for **ALL** columns. Default value is TRUE.
#
# Returns a data.frame or tibble, based on input table type, with additonal
# annotation columns
annotate_long_format_table <- function(long_format_table,
                                       is_gene_level_table = FALSE,
                                       add_Protein_RefSeq_ID = FALSE,
                                       replace_na_with_empty_string = TRUE) {
  # Check input long_format_table class is tibble or data.frame
  if (tibble::is_tibble(long_format_table)) {
    output_tbl_type <- "tibble"
  } else if (is.data.frame(long_format_table)) {
    output_tbl_type <- "data.frame"
  } else {
    stop(paste0(
      "Unsupported table type ", class(long_format_table),
      ".\nTry to convert the long_format_table to",
      " tibble::tibble or data.frame."))
  }
  # Check colnames(long_format_table) is character
  if (is.null(colnames(long_format_table))) {
    stop("colnames(long_format_table) is NULL. Check data integrity.")
  } else if (!is.character(colnames(long_format_table))) {
    stop("colnames(long_format_table) is not character. Check data integrity.")
  } 
  # Check required columns are in long_format_table walk is like sapply, but it
  # returns input argument invisibly even if there is any other return
  # values
  # Ref: https://purrr.tidyverse.org/reference/index.html
  purrr::walk(
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease"),
    function(req_col) {
      if (!(req_col %in% colnames(long_format_table))) {
        stop(paste0(
          req_col, " is not in colnames(long_format_table).\n",
          "Add a ", req_col, " column or rename an existing column to ",
          req_col, "."))
      }
    }
  )

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

  # Read annotation data and remove duplicated rows
  # Gene ENSG ID -> gene full name and protein refseq ID
  ensg_gname_prt_refseq_df <- dplyr::distinct(readr::read_tsv(
    file.path(root_dir, "analyses", "long-format-table-utils", "annotator",
              "annotation-data", "ensg-gene-full-name-refseq-protein.tsv"),
    col_types = readr::cols(.default = readr::col_guess()),
    guess_max = 100000))
  # Gene hugo symbol -> OncoKB cancer gene and oncogene/TSG
  hgsb_oncokb_cgene_oncogene_tsg_df <- dplyr::distinct(readr::read_tsv(
    file.path(root_dir, "analyses", "long-format-table-utils", "annotator",
              "annotation-data", "oncokb-cancer-gene-list.tsv"),
    col_types = readr::cols(.default = readr::col_guess()),
    guess_max = 100000))
  # Gene hugo symbol -> gene type
  hgsb_gtype_df <- dplyr::distinct(readr::read_tsv(
    file.path(root_dir, "analyses", "fusion_filtering", "references",
              "genelistreference.txt"),
    col_types = readr::cols(.default = readr::col_guess()),
    guess_max = 100000))
  # Gene ENSG ID -> gene hugo symbol, RMTL, and RMTL version
  ensg_hgsb_rmtl_df <- dplyr::distinct(readr::read_tsv(
    file.path(root_dir, "data", "ensg-hugo-rmtl-v1-mapping.tsv"),
    col_types = readr::cols(.default = readr::col_guess()),
    guess_max = 100000))
  # cancer_group -> EFO and MONDO
  cgroup_efo_mondo_df <- dplyr::distinct(readr::read_tsv(
    file.path(root_dir, "data", "efo-mondo-map.tsv"),
    col_types = readr::cols(.default = readr::col_guess()),
    guess_max = 100000))



}
