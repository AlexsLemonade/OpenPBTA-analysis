# Add gene and cancer_group annotations to a long-format table
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


  # Check no NA or duplicate in the key columns of annotation data tables
  #
  # For an annotation table, NA or duplicate should not exist in the key
  # column that is used to join_by for adding annotations to the input table

  # assert no NA in gene symbols
  if (!identical(sum(is.na(dplyr::pull(ensg_gname_prt_refseq_df,
                                       Gene_Ensembl_ID))),
                 as.integer(0))) {
    stop(paste0("analyses/long-format-table-utils/annotator/annotation-data/",
                "ensg-gene-full-name-refseq-protein.tsv ",
                "has NAs in the Gene_Ensembl_ID column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # assert all symbols are unique
  if (!identical(nrow(ensg_gname_prt_refseq_df),
                 length(unique(dplyr::pull(ensg_gname_prt_refseq_df,
                                           Gene_Ensembl_ID))))) {
    stop(paste0("analyses/long-format-table-utils/annotator/annotation-data/",
                "ensg-gene-full-name-refseq-protein.tsv ",
                "has duplicates in the Gene_Ensembl_ID column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }

  # assert no NA in gene symbols
  if (!identical(sum(is.na(dplyr::pull(hgsb_oncokb_cgene_oncogene_tsg_df,
                                       `Hugo Symbol`))),
                 as.integer(0))) {
    stop(paste0("analyses/long-format-table-utils/annotator/annotation-data/",
                "oncokb-cancer-gene-list.tsv ",
                "has NAs in the `Hugo Symbol` column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # assert all symbols are unique
  if (!identical(nrow(hgsb_oncokb_cgene_oncogene_tsg_df),
                 length(unique(dplyr::pull(hgsb_oncokb_cgene_oncogene_tsg_df,
                                           `Hugo Symbol`))))) {
    stop(paste0("analyses/long-format-table-utils/annotator/annotation-data/",
                "oncokb-cancer-gene-list.tsv ",
                "has duplicates in the `Hugo Symbol` column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }

  # assert no NA in hgsb_gtype_df
  if (!identical(sum(is.na(hgsb_gtype_df)), as.integer(0))) {
    stop(paste0("analyses/fusion_filtering/references/genelistreference.txt ",
                "has NAs in the table.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }

  # assert all ensg_ids and gene_symbols are not NA
  if (!identical(sum(is.na(ensg_hgsb_rmtl_df$ensg_id)), as.integer(0))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has NAs in the ensg_id column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  if (!identical(sum(is.na(ensg_hgsb_rmtl_df$gene_symbol)), as.integer(0))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has NAs in the gene_symbol column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # assert all ensg_id are unique
  if (!identical(length(unique(ensg_hgsb_rmtl_df$ensg_id)),
                 nrow(ensg_hgsb_rmtl_df))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has duplicates in the ensg_id column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # asert all rmtl NAs have version NAs, vice versa
  if (!identical(is.na(ensg_hugo_rmtl_df$rmtl),
                 is.na(ensg_hugo_rmtl_df$version))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has rmtl column NAs with version column non-NAs, or",
                "version column NAs with rmtl column non-NAs\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }

  # assert all cancer_groups are not NA
  if (!identical(sum(is.na(cgroup_efo_mondo_df$cancer_group)), as.integer(0))) {
    stop(paste0("efo-mondo-map.tsv has NAs in the cancer_group column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # assert all cancer_groups are unique.
  if (!identical(length(unique(cgroup_efo_mondo_df$cancer_group)),
                 nrow(cgroup_efo_mondo_df))) {
    stop(
      paste0("efo-mondo-map.tsv has duplicates in the cancer_group column.\n",
             "Check data integrity. Submit a data question GitHub issue."))
  }


}
