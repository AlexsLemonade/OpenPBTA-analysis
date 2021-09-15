# Add gene and cancer_group annotations to a long-format table
#
# Args:
# - long_format_table: A tibble that has zero or more of the following character
#   columns that are required for adding their corresponding annotation columns:
#   - Gene_symbol: HUGO symbols, e.g. PHLPP1, TM6SF1, and DNAH5. The Gene_symbol
#     column is required for adding the following annotation columns:
#     - Gene_type
#     - OncoKB_cancer_gene
#     - OncoKB_oncogene_TSG
#   - Gene_Ensembl_ID: Ensembl ENSG IDs without `.#` versions, e.g.
#     ENSG00000039139, ENSG00000111261, and ENSG00000169710. The Gene_Ensembl_ID
#     column is required for adding the following annotation columns:
#     - RMTL
#     - Gene_full_name
#     - Protein_RefSeq_ID
#   - Disease: The `cancer_group` in the `histologies.tsv`, e.g.
#     Adamantinomatous Craniopharyngioma, Atypical Teratoid Rhabdoid Tumor, and
#     Low-grade glioma/astrocytoma. The Disease column is required for adding
#     the following annotation columns:
#     - EFO
#     - MONDO
#   - GTEx_tissue_group: The `gtex_group` in the `histologies.tsv`, e.g.
#     Adipose, Kidney, and Thyroid. The GTEx_tissue_group column is required for
#     adding the following annotation column:
#     - GTEx_tissue_group_UBERON
#   - GTEx_tissue_subgroup: The `gtex_subgroup` in the `histologies.tsv`, e.g.
#     Adipose - Subcutaneous, Kidney - Cortex, and Thyroid. The
#     GTEx_tissue_subgroup column is required for adding the following
#     annotation column:
#     - GTEx_tissue_subgroup_UBERON
# - columns_to_add: a character vector of unique names of the columns to be
#   added to the input table. The vector can contain zero or more of the
#   following column names: "RMTL", "Gene_type", "OncoKB_cancer_gene",
#   "OncoKB_oncogene_TSG", "Gene_full_name", "Protein_RefSeq_ID", "EFO",
#   "MONDO", "GTEx_tissue_group_UBERON", "GTEx_tissue_subgroup_UBERON". Default
#   value is to add "RMTL", "Gene_type", "OncoKB_cancer_gene",
#   "OncoKB_oncogene_TSG", "Gene_full_name", "Protein_RefSeq_ID", "EFO", and
#   "MONDO".
#   - Notes:
#     - match.arg is **NOT** used for this parameter, because exact matches are
#       required
#     - The order of returned annotated table has the same column order as the
#       input table and added columns in the order of columns_to_add parameter
# - replace_na_with_empty_string: TRUE or FALSE on whether NAs should be
#   replaced with empty strings for **ALL COLUMNS THAT HAVE NA** in the output
#   table. Default value is TRUE.
#
# Returns a tibble with additonal annotation columns
#
# Notes on requiring Gene_symbol and Gene_Ensembl_ID:
# - Some Gene_symbols are mapped to multiple Gene_Ensembl_IDs, so adding
#   Gene_Ensembl_IDs by mapping Gene_symbols with
#   data/ensg-hugo-rmtl-mapping.tsv may implicitly introduce duplicated rows.
#   Therefore, adding Gene_Ensembl_IDs by mapping Gene_symbols is left to users
#   with cautions for potentially introducing unwanted duplicates.
# - Similarly, some Gene_Ensembl_IDs are mapped to multiple Gene_symbols, so
#   adding Gene_symbols by mapping Gene_Ensembl_IDs with
#   data/ensg-hugo-rmtl-mapping.tsv may implicitly introduce duplicated rows.
#   Therefore, adding Gene_symbols by mapping Gene_Ensembl_IDs is left to users
#   with cautions for potentially introducing unwanted duplicates.
# - Certain annotation files use Gene_symbol as key columns, and certain other
#   annotation files use Gene_Ensembl_ID as key columns
#
# Notes on how to add new annotation columns for maintainers:
# - Update this documentation comment for the long_format_table and
#   columns_to_add parameters
# - For backward compatibility, do not change the default parameter value of
#   columns_to_add
# - Add new column names to the available_ann_columns variable
# - Implement annotation procedures to add new columns by adding additional
#   ann_tbl_l elements
# - Update annotator-cli.R help message for -c/--columns-to-add
# - Add new tests to tests/test_annotate_long_format_table.R and
#   tests/test_annotator_cli.R for new annotation columns
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
                             "Protein_RefSeq_ID", "EFO", "MONDO",
                             "GTEx_tissue_group_UBERON",
                             "GTEx_tissue_subgroup_UBERON")
  # Helper function to convert a character vector to a single string
  char_vec_to_str <- function(x) {
    return(paste0("c(\"", paste(x, collapse = "\", \""), "\")"))
  }
  # Assert columns_to_add columns are available and not already in the input
  # long_format_table
  #
  # walk is like sapply, but it returns input argument invisibly even if there
  # is any other return values
  #
  # Ref:
  # https://purrr.tidyverse.org/reference/index.html
  purrr::walk(
    columns_to_add,
    function(col_add) {
      if (!(col_add %in% available_ann_columns)) {
        stop(paste0(
          col_add, " is not available.\n",
          "Available annotation columns are ",
          char_vec_to_str(available_ann_columns), "."))
      }
      if (col_add %in% colnames(long_format_table)) {
        stop(paste0(col_add, " is already in the input long_format_table."))
      }
    }
  )
  # Keep the order of columns for the returned table
  ret_tbl_col_order <- c(colnames(long_format_table), columns_to_add)

  # Detect the ".git" folder -- this will be in the project root directory. Use
  # this as the root directory to ensure proper execution, no matter where it is
  # called from.
  #
  # This only works if the working directory is OpenPedCan-analysis or a
  # subdirectory of OpenPedCan-analysis
  #
  # root_dir is the absolute path of OpenPedCan-analysis
  #
  # Adapted from the oncoprint-landscape module
  #
  # rprojroot::has_file(".git/index") returns a rprojroot::root_criterion, and
  # main git working tree, created by git clone and git init, has the .git/index
  # file
  #
  # rprojroot::has_file(".git") returns a rprojroot::root_criterion, and linked
  # git working tree, created by git worktree add, has the .git file
  #
  # "Root criteria can be combined with the | operator. The result is a
  # composite root criterion that requires either of the original criteria to
  # match." -- help("root_criterion", "rprojroot") rprojroot_1.3-2
  tryCatch(
    {
      root_dir <- rprojroot::find_root(
        rprojroot::has_file(".git/index") | rprojroot::has_file(".git"))
    },
    error = function(err_cond) {
      # adapted from http://adv-r.had.co.nz/Exceptions-Debugging.html
      err_cond$message <- paste0(
        err_cond$message,
        "\nTry re-running this function with working directory as ",
        "OpenPedCan-analysis or a subdirectory of OpenPedCan-analysis.\n")
      stop(err_cond)
    }
  )

  # List of required annotation table file paths, join_by columns, and tables
  #
  # Each element of ann_tbl_l is a list(file_path = file_path, join_by_column =
  # join_by_column, tibble = tibble). The tibble value is the table that is read
  # in, after the list of paths is created. The tibble value may be further
  # processed before annotating the input long_format_table.
  #
  # NOTE: join_by_column must be a length 1 character vector that has no name
  # attribute, which is passed to left_join(long_format_table, annotation_table,
  # by = join_by_column). If the annotation table join_by_column has a different
  # name, rename it to be the same as the one in long_format_table.
  #
  # NA or duplicate should not exist in the join_by_column.
  #
  # Steps:
  # - Initialize a list that contains annotation table path, join_by_column, and
  #   the annotation tibble that is read in. The join_by_column have to be
  #   there, so proper checks can be performed in a functional way.
  # - Process annotation table, and the processed table should have the
  #   join_by_column. Again, the join_by_column should exist in both
  #   long_format_table and annotation_table tables.

  # annotation table list
  ann_tbl_l <- list()

  # Helper function to read all columns as characters, as no
  # non-character/string operation is needed, and remove duplicated rows
  #
  # Args:
  # - file_path: a sinigle character path to the annotation table
  #
  # Returns a tibble that is read in from the path
  rmrdup_all_char_col_read_tsv <- function(file_path) {
    ann_tibble <- dplyr::distinct(
      readr::read_tsv(
        file_path,
        col_types = readr::cols(.default = readr::col_character())
      )
    )
    return(ann_tibble)
  }
  # Helper function to check the type of join_by_column value in ann_tbl_l
  # element
  #
  # Args:
  # - join_by_column: a singele character value of the column to be joined by
  #
  # Returns NULL, as the function only checks whether the input join_by_column
  # is valid or not
  check_join_by_column <- function(join_by_column) {
    # Check join_by_column is valid
    if (is.null(join_by_column)) {
      stop(paste0("annotate_long_format_table function internal error. ",
                  "Submit a GitHub Issue with your error ",
                  "message and input data.\n",
                  "join_by_column cannot be NULL."))
    }
    if (!(is.character(join_by_column) &&
          identical(length(join_by_column), as.integer(1)) &&
          is.null(names(join_by_column)))) {
      stop(paste0("annotate_long_format_table function internal error. ",
                  "Submit a GitHub Issue with the full error traceback",
                  "message and input data.\n",
                  "join_by_column for must be a single character value ",
                  "that has no name attribute."))
    }
    return(NULL)
  }

  # Helper function to initialize ann_tbl_l element, by adding file_path,
  # join_by_column, and read table
  #
  # Args:
  # - file_path: a sinigle character path to the annotation table
  # - join_by_column: a singele character value of the column to be joined by
  #
  # Returns a list(file_path = file_path, join_by_column = join_by_column,
  # tibble = tibble_that_is_read_in)
  init_ann_tbl_l_element <- function(file_path, join_by_column) {
    check_join_by_column(join_by_column)

    ann_tbl_l_element <- list(
      file_path = file_path, join_by_column = join_by_column)

    ann_tbl_l_element$tibble <- rmrdup_all_char_col_read_tsv(
      ann_tbl_l_element$file_path)

    return(ann_tbl_l_element)
  }

  # Only add required annotation tables
  if (any(c("Protein_RefSeq_ID", "Gene_full_name") %in% columns_to_add)) {
    ann_tbl_l$mygene_info <- init_ann_tbl_l_element(
      file_path = file.path(
        root_dir, "analyses", "long-format-table-utils", "annotator",
        "annotation-data", "ensg-gene-full-name-refseq-protein.tsv"),
      join_by_column = "Gene_Ensembl_ID")
  }

  if (any(c("OncoKB_cancer_gene", "OncoKB_oncogene_TSG") %in% columns_to_add)) {
    ann_tbl_l$onkokb <- init_ann_tbl_l_element(
      file_path = file.path(
        root_dir, "analyses", "long-format-table-utils", "annotator",
        "annotation-data", "oncokb-cancer-gene-list.tsv"),
      join_by_column = "Gene_symbol")
    # Process tibble for joining
    ann_tbl_l$onkokb$tibble <- ann_tbl_l$onkokb$tibble %>%
      dplyr::select(
        `Hugo Symbol`, `Is Oncogene`, `Is Tumor Suppressor Gene`) %>%
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
  }

  if ("Gene_type" %in% columns_to_add) {
    ann_tbl_l$gene_type <- init_ann_tbl_l_element(
      file_path = file.path(
        root_dir, "analyses", "fusion_filtering", "references",
        "genelistreference.txt"),
      join_by_column = "Gene_symbol"
    )
    # Process tibble for joining
    ann_tbl_l$gene_type$tibble <- ann_tbl_l$gene_type$tibble %>%
      dplyr::group_by(Gene_Symbol) %>%
      dplyr::summarise(
        Gene_type = paste(sort(unique(purrr::discard(type, is.na))),
                          collapse = ","))  %>%
      dplyr::rename(Gene_symbol = Gene_Symbol)
  }

  if ("RMTL" %in% columns_to_add) {
    ann_tbl_l$fda_rmtl <- init_ann_tbl_l_element(
      file_path = file.path(root_dir, "data", "ensg-hugo-rmtl-mapping.tsv"),
      join_by_column = "Gene_Ensembl_ID")
    # Asert all rmtl NAs have version NAs, vice versa
    if (!identical(is.na(ann_tbl_l$fda_rmtl$tibble$rmtl),
                  is.na(ann_tbl_l$fda_rmtl$tibble$version))) {
      stop(paste0("ensg-hugo-rmtl-mapping.tsv ",
                  "has rmtl column NAs with version column non-NAs, or",
                  "version column NAs with rmtl column non-NAs\n",
                  "Check data integrity. Submit a data question GitHub issue."))
    }
    # Process tibble for joining
    ann_tbl_l$fda_rmtl$tibble <- ann_tbl_l$fda_rmtl$tibble %>%
      dplyr::select(ensg_id, rmtl, version) %>%
      dplyr::filter(!is.na(rmtl), !is.na(version)) %>%
      dplyr::mutate(RMTL = paste0(rmtl, " (", version, ")")) %>%
      dplyr::select(ensg_id, RMTL) %>%
      dplyr::rename(Gene_Ensembl_ID = ensg_id) %>%
      dplyr::distinct()
  }

  if (any(c("EFO", "MONDO") %in% columns_to_add)) {
    ann_tbl_l$cgroup_ontology = init_ann_tbl_l_element(
      file_path = file.path(root_dir, "data", "efo-mondo-map.tsv"),
      join_by_column = "Disease")
    # Process tibble for joining
    ann_tbl_l$cgroup_ontology$tibble <- dplyr::rename(
      ann_tbl_l$cgroup_ontology$tibble,
      Disease = cancer_group, EFO = efo_code, MONDO = mondo_code)
  }

  if ("GTEx_tissue_group_UBERON" %in% columns_to_add) {
    ann_tbl_l$gtex_group_ontology = init_ann_tbl_l_element(
      file_path = file.path(root_dir, "data", "uberon-map-gtex-group.tsv"),
      join_by_column = "GTEx_tissue_group")
    # Process tibble for joining
    ann_tbl_l$gtex_group_ontology$tibble <- dplyr::rename(
      ann_tbl_l$gtex_group_ontology$tibble,
      GTEx_tissue_group = gtex_group, GTEx_tissue_group_UBERON = uberon_code)
  }

  if ("GTEx_tissue_subgroup_UBERON" %in% columns_to_add) {
    ann_tbl_l$gtex_subgroup_ontology = init_ann_tbl_l_element(
      file_path = file.path(root_dir, "data", "uberon-map-gtex-subgroup.tsv"),
      join_by_column = "GTEx_tissue_subgroup")
    # Process tibble for joining
    ann_tbl_l$gtex_subgroup_ontology$tibble <- dplyr::rename(
      ann_tbl_l$gtex_subgroup_ontology$tibble,
      GTEx_tissue_subgroup = gtex_subgroup,
      GTEx_tissue_subgroup_UBERON = uberon_code)
  }

  # Check required columns are in long_format_table
  required_columns <- purrr::map_chr(ann_tbl_l, function(x) x$join_by_column)
  purrr::walk(
    required_columns,
    function(req_col) {
      if (!(req_col %in% colnames(long_format_table))) {
        stop(paste0(
          req_col, " is not in colnames(long_format_table).\n",
          "Add a ", req_col, " column or rename an existing column to ",
          req_col, ". Required columns are ",
          char_vec_to_str(required_columns), " to add ",
          char_vec_to_str(columns_to_add), " annotation columns."))
      }
      if (!is.character(long_format_table[[req_col]])) {
        stop(paste0(
          req_col, " is not a character column in the long_format_table.\n",
          "Convert the", req_col, " column to character."))
      }
    }
  )

  # Check no NA or duplicate in the join_by_column of annotation data tables
  purrr::walk(
    ann_tbl_l,
    function(xl) {
      # assert no NA in xl$join_by_column
      if (!identical(sum(is.na(xl$tibble[, xl$join_by_column])),
                     as.integer(0))) {
        stop(paste0(xl$file_path, "has NA(s) in the ",
                    xl$join_by_column, " column.\n",
                    "Check data integrity. Submit a data question ",
                    "GitHub issue."))
      }
      # assert no duplicate in xl$join_by_column
      if (!identical(nrow(xl$tibble),
                     nrow(dplyr::distinct(xl$tibble[, xl$join_by_column])))) {
        stop(paste0(xl$file_path, "has duplicate(s) in the ",
                    xl$join_by_column, " column.\n",
                    "Check data integrity. Submit a data question ",
                    "GitHub issue."))
      }
    }
  )


  # Annotate required columns for the input table
  #
  # .dir = "forward", so:
  # - x starts with .init value
  # - y is always the next element in the ann_tbl_l
  # - the function return value will be the next x
  ann_long_format_table <- purrr::reduce(
    ann_tbl_l,
    function(x, y) {
      # If one or more columns in y$tibble are already in x$tibble, remove them
      # before joining. The y$join_by_column should exist in both x$tibble and
      # y$tibble, so it should not be in y_cols_not_in_x, hence
      # rm_cmn_col_y_tibble will not have duplicated columns.
      #
      # The type of y$join_by_column is checked above, but check here again in
      # case it is modified by prior procedures
      check_join_by_column(y$join_by_column)

      y_cols_not_in_x <- colnames(y$tibble)[
        !colnames(y$tibble) %in% colnames(x$tibble)]

      if (!(is.character(colnames(x$tibble)) &&
            is.character(colnames(y$tibble)) &&
            (y$join_by_column %in% colnames(y$tibble)) &&
            (y$join_by_column %in% colnames(x$tibble)) &&
            (y_cols_not_in_x %in% colnames(y$tibble)) &&
            (!y_cols_not_in_x %in% colnames(x$tibble)) &&
            (!y$join_by_column %in% y_cols_not_in_x))) {
        stop(paste0("annotate_long_format_table function internal error. ",
                    "Submit a GitHub Issue with the full error traceback",
                    "message and input data.\n",
                    "join_by_column issue for ", y$file_path,
                    ". Inspect the joining part in ",
                    "annotate_long_format_table."))
      }
      rm_cmn_col_y_tibble <- y$tibble[, c(y$join_by_column, y_cols_not_in_x)]

      res_tibble <- dplyr::left_join(
        x$tibble, rm_cmn_col_y_tibble, by = y$join_by_column)

      return(list(tibble = res_tibble))
    },
    .init = list(tibble = long_format_table),
    .dir = "forward")$tibble

  # The OncoKB_cancer_gene values are Y or N only
  ann_long_format_table <- tidyr::replace_na(
    ann_long_format_table, list(OncoKB_cancer_gene = "N"))

  if (replace_na_with_empty_string) {
    # Only replace NA with empty string in columns that have NA in them. This
    # will not change the value types of the columns that have no NA, so the
    # output JSON/JSONL table will be backward compatible with the ones that
    # were generated without using the annotator API.
    ann_long_format_table <- dplyr::mutate_if(
      ann_long_format_table,
      function(x) sum(is.na(x)) > 0,
      function(x) tidyr::replace_na(x, replace = ""))
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
  # tidyselect::all_of, the docker image only has one_of
  ann_long_format_table <- dplyr::select(
    ann_long_format_table, tidyselect::one_of(ret_tbl_col_order))
  return(ann_long_format_table)
}
