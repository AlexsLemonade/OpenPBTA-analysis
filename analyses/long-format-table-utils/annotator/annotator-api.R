# Add gene and cancer_group annotations to a long-format table
#
# Args:
# - long_format_table: A tibble that contains the following columns as
#   character:
#   - Gene_symbol: HUGO symbols, e.g. PHLPP1, TM6SF1, and DNAH5
#   - Gene_Ensembl_ID: Ensembl ENSG IDs without `.#` versions, e.g.
#     ENSG00000039139, ENSG00000111261, and ENSG00000169710
#   - Disease: The `cancer_group` in the `histologies.tsv`, e.g.
#     Adamantinomatous Craniopharyngioma, Atypical Teratoid Rhabdoid Tumor, and
#     Low-grade glioma/astrocytoma
# - columns_to_add: a character vector of unique names of the columns to be
#   added to the input table. The vector can contain zero or more of the
#   following column names: "RMTL", "Gene_type", "OncoKB_cancer_gene",
#   "OncoKB_oncogene_TSG", "Gene_full_name", "Protein_RefSeq_ID", "EFO", and
#   "MONDO". Default value is to add all columns. Notes:
#   - match.arg is **NOT** used for this parameter, because exact matches are
#     required
#   - The order of returned annotated table has the same column order as the
#     input table and added columns in the order of columns_to_add parameter.
# - replace_na_with_empty_string: TRUE or FALSE on whether to replace NAs with
#   empty strings for **ALL** columns of the input table. Default value is TRUE.
#
# Returns a tibble with additonal annotation columns
annotate_long_format_table <- function(
  long_format_table,
  columns_to_add = c("RMTL", "Gene_type", "OncoKB_cancer_gene",
                     "OncoKB_oncogene_TSG", "Gene_full_name",
                     "Protein_RefSeq_ID", "EFO", "MONDO"),
  replace_na_with_empty_string = TRUE) {
  # Get %>% without loading the whole library
  `%>%` <- magrittr::`%>%`

  # Check input long_format_table class is tibble
  #
  # The function now only supports tibble, because handling data.frame or other
  # types of table is complex.
  #
  # Support for other types of table could be added if required at a later
  # point.
  if (!tibble::is_tibble(long_format_table)) {
    stop(paste0(
      "Unsupported table type ", class(long_format_table),
      ".\nTry converting the long_format_table to",
      " tibble::tibble with tibble::as_tibble."))
  }
  # Check colnames(long_format_table) is character
  if (is.null(colnames(long_format_table))) {
    stop("colnames(long_format_table) is NULL. Check data integrity.")
  } else if (!is.character(colnames(long_format_table))) {
    stop("colnames(long_format_table) is not character. Check data integrity.")
  }
  # Helper function to convert a character vector to a single string
  char_vec_to_str <- function(x) {
    return(paste0("c(\"", paste(x, collapse = "\", \""), "\")"))
  }
  # Check required columns are in long_format_table walk is like sapply, but it
  # returns input argument invisibly even if there is any other return
  # values
  # Ref: https://purrr.tidyverse.org/reference/index.html
  required_columns <- c("Gene_symbol", "Gene_Ensembl_ID", "Disease")
  purrr::walk(
    required_columns,
    function(req_col) {
      if (!(req_col %in% colnames(long_format_table))) {
        stop(paste0(
          req_col, " is not in colnames(long_format_table).\n",
          "Add a ", req_col, " column or rename an existing column to ",
          req_col, ". Required columns are ",
          char_vec_to_str(required_columns), "."))
      }
    }
  )
  # The columns_to_add must match exactly, so match.arg is not used.
  if (is.null(columns_to_add)) {
    stop(paste0("columns_to_add cannot be NULL. ",
                "Read the documentation/comment/docstring for usage."))
  }
  if (!is.character(columns_to_add)) {
    stop(paste0("columns_to_add must be a character vector. ",
                "Read the documentation/comment/docstring for usage."))
  }
  if (!identical(length(columns_to_add), length(unique(columns_to_add)))) {
    stop(paste0("columns_to_add values must be unique. ",
                "Read the documentation/comment/docstring for usage."))
  }
  if (identical(columns_to_add, character(0))) {
    # no column to add, so return the input table
    return(long_format_table)
  }
  available_ann_columns <- c("RMTL", "Gene_type", "OncoKB_cancer_gene",
                             "OncoKB_oncogene_TSG", "Gene_full_name",
                             "Protein_RefSeq_ID", "EFO", "MONDO")
  purrr::walk(
    columns_to_add,
    function(col_add) {
      if (!(col_add %in% available_ann_columns)) {
        stop(paste0(
          col_add, " is not available.\n",
          "Available annotation columns are ",
          char_vec_to_str(available_ann_columns), "."))
      }
    }
  )
  # Keep the order of columns for the returned table
  ret_tbl_col_order <- c(colnames(long_format_table), columns_to_add)

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
  # Gene ENSG ID -> RMTL
  ensg_rmtl_df <- dplyr::distinct(readr::read_tsv(
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
  if (!identical(sum(is.na(ensg_rmtl_df$ensg_id)), as.integer(0))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has NAs in the ensg_id column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # assert all ensg_id are unique
  if (!identical(length(unique(ensg_rmtl_df$ensg_id)), nrow(ensg_rmtl_df))) {
    stop(paste0("ensg-hugo-rmtl-v1-mapping.tsv ",
                "has duplicates in the ensg_id column.\n",
                "Check data integrity. Submit a data question GitHub issue."))
  }
  # asert all rmtl NAs have version NAs, vice versa
  if (!identical(is.na(ensg_rmtl_df$rmtl),
                 is.na(ensg_rmtl_df$version))) {
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


  # Process annotation data for joining input table
  #
  # pp in variable names means preprocessed
  #
  # %>% is not used, so the user's environment will not be modified by importing
  # the %>% in this file
  pp_hgsb_oncokb_cgene_oncogene_tsg_df <- hgsb_oncokb_cgene_oncogene_tsg_df %>%
    dplyr::select(`Hugo Symbol`, `Is Oncogene`, `Is Tumor Suppressor Gene`) %>%
    dplyr::rename(
      Gene_symbol = `Hugo Symbol`, is_onco = `Is Oncogene`,
      is_tsg = `Is Tumor Suppressor Gene`) %>%
    dplyr::mutate(
      OncoKB_cancer_gene = "Y",
      OncoKB_oncogene_TSG = dplyr::case_when(
        is_onco == "Yes" & is_tsg == "Yes" ~ "Oncogene,TumorSuppressorGene",
        is_onco == "Yes" ~ "Oncogene",
        is_tsg == "Yes" ~ "TumorSuppressorGene",
        TRUE ~ as.character(NA))) %>%
    dplyr::select(Gene_symbol, OncoKB_cancer_gene, OncoKB_oncogene_TSG)

  pp_hgsb_gtype_df <- hgsb_gtype_df %>%
    dplyr::group_by(Gene_Symbol) %>%
    dplyr::summarise(
      Gene_type = paste(sort(unique(purrr::discard(type, is.na))),
                        collapse = ",")) %>%
    dplyr::rename(Gene_symbol = Gene_Symbol)

  pp_ensg_rmtl_df <- ensg_rmtl_df %>%
    dplyr::select(ensg_id, rmtl, version) %>%
    dplyr::filter(!is.na(rmtl), !is.na(version)) %>%
    dplyr::mutate(RMTL = paste0(rmtl, " (", version, ")")) %>%
    dplyr::select(ensg_id, RMTL) %>%
    dplyr::rename(Gene_Ensembl_ID = ensg_id)

  pp_cgroup_efo_mondo_df <- dplyr::rename(
    cgroup_efo_mondo_df,
    Disease = cancer_group, EFO = efo_code, MONDO = mondo_code)

  # # Code for checking data
  # print(head(ensg_gname_prt_refseq_df))
  # print(head(pp_hgsb_oncokb_cgene_oncogene_tsg_df))
  # print(head(pp_hgsb_gtype_df))
  # print(head(pp_ensg_rmtl_df))
  # print(head(pp_cgroup_efo_mondo_df))


  # Annotate required columns for the input table
  ann_long_format_table <- long_format_table
  if (any(c("Protein_RefSeq_ID", "Gene_full_name") %in% columns_to_add)) {
    ann_long_format_table <- dplyr::left_join(
      ann_long_format_table, ensg_gname_prt_refseq_df, by = "Gene_Ensembl_ID")
  }

  if (any(c("OncoKB_cancer_gene", "OncoKB_oncogene_TSG") %in% columns_to_add)) {
    ann_long_format_table <- dplyr::left_join(
      ann_long_format_table, pp_hgsb_oncokb_cgene_oncogene_tsg_df,
      by = "Gene_symbol")
    ann_long_format_table <- tidyr::replace_na(
      ann_long_format_table, list(OncoKB_cancer_gene = "N"))
  }

  if ("Gene_type" %in% columns_to_add) {
    ann_long_format_table <- dplyr::left_join(
      ann_long_format_table, pp_hgsb_gtype_df, by = "Gene_symbol")
  }

  if ("RMTL" %in% columns_to_add) {
    ann_long_format_table <- dplyr::left_join(
      ann_long_format_table, pp_ensg_rmtl_df, by = "Gene_Ensembl_ID")
  }

  if (any(c("EFO", "MONDO") %in% columns_to_add)) {
    ann_long_format_table <- dplyr::left_join(
      ann_long_format_table, pp_cgroup_efo_mondo_df, by = "Disease")
  }

  if (replace_na_with_empty_string) {
    ann_long_format_table <- dplyr::mutate_all(
      ann_long_format_table, function(x) tidyr::replace_na(x, replace = ""))
  }

  # Reorder the columns of the returned table to have the same column order as
  # the input table and added columns in the order of columns_to_add parameter
  #
  # Assert that input table columns and columns_to_add columns are in the
  # ann_long_format_table
  if (!all(ret_tbl_col_order %in% colnames(ann_long_format_table))) {
    stop(paste0("annotate_long_format_table function internal error. ",
                "Submit a GitHub Issue with your error ",
                "message and input data.\n",
                "ann_long_format_table column names are ",
                char_vec_to_str(colnames(ann_long_format_table)),
                ", which do not include all column names to be returned ",
                char_vec_to_str(ret_tbl_col_order), "."))
  }
  # Though tidyselect::one_of is superseded in favor of tidyselect::any_of and
  # tidyselect::all_of, the docker image only has one_of.
  ann_long_format_table <- dplyr::select(
    ann_long_format_table, tidyselect::one_of(ret_tbl_col_order))
  return(ann_long_format_table)
}
