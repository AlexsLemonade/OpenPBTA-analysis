# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_annotator_cli.R")


working_input_tsv_path <- "test_data/test_long_format_table.tsv"

# Add [] after reading to be compatible with readr >= 1.3.1, otherwise the tests
# will fail on readr >= 1.3.1 as found by @NHJohnson at
# <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/56
#  #issuecomment-885188592>
#
# readr 1.3.1 returns spec_tbl_df subclass, which becomes tbl_df after any
# subsetting
#
# Ref: https://www.tidyverse.org/blog/2018/12/readr-1-3-1/#tibble-subclass
long_format_tibble <- readr::read_tsv(
  working_input_tsv_path,
  col_types = readr::cols(.default = readr::col_character()))[]

inspected_annotated_long_format_tibble <- readr::read_tsv(
  "test_data/inspected_annotated_test_long_format_table.tsv",
  col_types = readr::cols(.default = readr::col_character()),
  na = c("NA"), quoted_na = FALSE, trim_ws = FALSE)[]

# esana means empty string as NA. This table is used to test the -r option
inspected_annotated_long_format_tibble_esana <- readr::read_tsv(
  "test_data/inspected_annotated_test_long_format_table.tsv",
  col_types = readr::cols(.default = readr::col_character()),
  na = c("", "NA"), quoted_na = FALSE, trim_ws = FALSE)[]

testthat::expect_gt(
  sum(is.na(inspected_annotated_long_format_tibble_esana)), 0)

testthat::expect_equal(
  sum(is.na(inspected_annotated_long_format_tibble)), 0)

testthat::expect_equal(
  sum(is.na(long_format_tibble)), 0)

# to save intermediate files
scratch_dir <- "test_scratch"
annotator_cli_path <- file.path("..", "annotator-cli.R")
annotator_cli_output_path <- file.path(
  scratch_dir, "annotated_test_long_format_table.tsv")

# Test cases
#
# Helper function to run annotator CLI
run_cli_get_tibble <- function(columns_to_add,
                               input_table_path,
                               output_table_path,
                               remove_input_table = FALSE,
                               replace_na_with_empty_string = TRUE,
                               read_output_as_a_string = FALSE) {
  run_command <- function(cmd_str) {
    # ignore.stderr = TRUE avoids printing when testing
    #
    # intern = TRUE captures output as return value
    #
    # If a command fails, only a warning is generated, so expect warning when
    # testing
    #
    # Future updates could parse warning messages to determine whether it is
    # failed by the CLI or other issues
    out <- system(cmd_str, ignore.stderr = TRUE, intern = TRUE)
    return(out)
  }
  cmd_pfx <- "Rscript --vanilla ../annotator-cli.R "
  if (replace_na_with_empty_string) {
    cmd_pfx <- paste0(cmd_pfx, "-r ")
  }

  if (is.null(columns_to_add)) {
    # run without -c
    run_command(paste0(
      cmd_pfx,
      " -i ", input_table_path,
      " -o ", output_table_path))
  } else {
    columns_to_add_opt_val <- paste0(
      "'", paste(columns_to_add, collapse = ","), "'")
    run_command(paste0(
      cmd_pfx,
      "-c ", columns_to_add_opt_val,
      " -i ", input_table_path,
      " -o ", output_table_path))
  }

  if (remove_input_table) {
    # Input table may not exist, which is to test fail
    #
    # clean up, so other tests will not be affected
    if (file.exists(input_table_path)) {
      file.remove(input_table_path)
    }
  }

  # the file may not be created due to CLI call failure
  if (file.exists(output_table_path)) {
    if (read_output_as_a_string) {
      cli_res <- readr::read_file(output_table_path)
    } else {
      cli_res <- readr::read_tsv(
        output_table_path,
        col_types = readr::cols(.default = readr::col_character()),
        na = c("NA"), quoted_na = FALSE, trim_ws = FALSE)[]
    }

    # clean up, so other tests will not be affected
    file.remove(output_table_path)
    return(cli_res)
  } else {
    return(NULL)
  }
}

# Add package prefix for auto completion purpose only
#
# Test standard usecase
testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble)

# The -r/--replace-na-with-empty-string does not need to be tested for column
# type conversions. The CLI write_tsv output table is the same even if a
# non-character column is converted to a character column before write_tsv,
# because "Values are only quoted if they contain a comma, quote or newline" (--
# help("write_tsv", "readr") 1.3.1).
#
# run_cli_get_tibble only reads "" as "" for character columns
testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path,
    replace_na_with_empty_string = FALSE),
  inspected_annotated_long_format_tibble_esana)

testthat::expect_true(
  stringr::str_detect(
    run_cli_get_tibble(
      columns_to_add = NULL,
      input_table_path = working_input_tsv_path,
      output_table_path = annotator_cli_output_path,
      replace_na_with_empty_string = FALSE,
      read_output_as_a_string = TRUE),
    "\tNA\t"))

testthat::expect_false(
  stringr::str_detect(
    run_cli_get_tibble(
      columns_to_add = NULL,
      input_table_path = working_input_tsv_path,
      output_table_path = annotator_cli_output_path,
      replace_na_with_empty_string = TRUE,
      read_output_as_a_string = TRUE),
    "\tNA\t"))

testthat::expect_true(
  stringr::str_detect(
    run_cli_get_tibble(
      columns_to_add = NULL,
      input_table_path = working_input_tsv_path,
      output_table_path = annotator_cli_output_path,
      replace_na_with_empty_string = TRUE,
      read_output_as_a_string = TRUE),
    "\t\t"))

testthat::expect_false(
  stringr::str_detect(
    run_cli_get_tibble(
      columns_to_add = NULL,
      input_table_path = working_input_tsv_path,
      output_table_path = annotator_cli_output_path,
      replace_na_with_empty_string = FALSE,
      read_output_as_a_string = TRUE),
    "\t\t"))

# Test annotation order
testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("RMTL", "Gene_type", "OncoKB_cancer_gene",
                       "OncoKB_oncogene_TSG", "Gene_full_name",
                       "Protein_RefSeq_ID", "EFO", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean", "RMTL",
      "Gene_type", "OncoKB_cancer_gene", "OncoKB_oncogene_TSG",
      "Gene_full_name", "Protein_RefSeq_ID", "EFO", "MONDO")])

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("MONDO", "RMTL", "EFO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "MONDO", "RMTL", "EFO")])

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("RMTL", "EFO", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "RMTL", "EFO", "MONDO")])

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("RMTL", "Protein_RefSeq_ID", "Gene_full_name"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "RMTL", "Protein_RefSeq_ID", "Gene_full_name")])

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("OncoKB_oncogene_TSG", "OncoKB_cancer_gene", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "OncoKB_oncogene_TSG", "OncoKB_cancer_gene", "MONDO")])

# Return same table if no annotation to add
testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = "",
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path),
  long_format_tibble)

# Error on duplicated annotation columns
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("OncoKB_oncogene_TSG", "OncoKB_oncogene_TSG", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("OncoKB_oncogene_TSG", "MONDO", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

# Error on non-available annotation columns
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("NOT_AVAILABLE", "MONDO", "MONDO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("NOT_AVAILABLE"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

# Error on missing required columns
req_col_missing_tbl_path <- file.path(
  scratch_dir, "test_missing_req_col_long_format_table.tsv")

readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_symbol),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(long_format_tibble, -Disease),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_symbol, -Gene_Ensembl_ID),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_symbol, -Gene_Ensembl_ID, -Disease),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID, -Disease),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_symbol),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("Gene_type"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Disease),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("EFO"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("Gene_full_name"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID),
  req_col_missing_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("RMTL", "Gene_full_name"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

# No error if all required columns are provided
readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID, -Disease),
  req_col_missing_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("Gene_type"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -Gene_Ensembl_ID, -Disease,
    -OncoKB_cancer_gene, -OncoKB_oncogene_TSG, -RMTL, -Gene_full_name,
    -Protein_RefSeq_ID, -EFO, -MONDO))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Gene_Ensembl_ID, -Gene_symbol),
  req_col_missing_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("EFO"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -Gene_Ensembl_ID, -Gene_symbol,
    -Gene_type, -OncoKB_cancer_gene, -OncoKB_oncogene_TSG, -RMTL,
    -Gene_full_name, -Protein_RefSeq_ID, -MONDO))

readr::write_tsv(
  dplyr::select(long_format_tibble, -Disease, -Gene_symbol),
  req_col_missing_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("Gene_full_name", "RMTL"),
    input_table_path = req_col_missing_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -Disease, -Gene_symbol,
    -Gene_type, -OncoKB_cancer_gene, -OncoKB_oncogene_TSG,
    -Protein_RefSeq_ID, -EFO, -MONDO,
    -Gene_full_name, -RMTL,
    Gene_full_name, RMTL))

# Error on requiring existing annotation columns
ann_col_exist_tbl_path <- file.path(
  scratch_dir, "test_ann_col_exist_long_format_table.tsv")

readr::write_tsv(
  dplyr::select(inspected_annotated_long_format_tibble, -EFO),
  ann_col_exist_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = ann_col_exist_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))

readr::write_tsv(
  inspected_annotated_long_format_tibble,
  ann_col_exist_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = NULL,
    input_table_path = ann_col_exist_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


readr::write_tsv(
  dplyr::select(inspected_annotated_long_format_tibble, -EFO),
  ann_col_exist_tbl_path)

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "MONDO"),
    input_table_path = ann_col_exist_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE))


# Error on duplicated annotation columns
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "EFO"),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

# Error on non character annotation columns
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c(1),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))

testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c(TRUE),
    input_table_path = working_input_tsv_path,
    output_table_path = annotator_cli_output_path))


# No error on requiring non-existing annotation columns
#
# Relocate to last, so the order is expected. Adapted from
# https://stackoverflow.com/a/43902237/4638182. The dplyr::relocate is not
# available in the Docker image
#
# The behavior of testthat::expect_equal changed at some point. The Docker
# image/container version does not check column order, whereas the latest
# version checks.
req_non_existing_ann_tbl_path <- file.path(
  scratch_dir, "test_req_non_existing_ann_long_format_table.tsv")

readr::write_tsv(
  dplyr::select(
    inspected_annotated_long_format_tibble, -EFO, -OncoKB_cancer_gene),
  req_non_existing_ann_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "OncoKB_cancer_gene"),
    input_table_path = req_non_existing_ann_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -EFO, -OncoKB_cancer_gene,
    EFO, OncoKB_cancer_gene))


readr::write_tsv(
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -EFO, -OncoKB_cancer_gene, -OncoKB_oncogene_TSG),
  req_non_existing_ann_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "OncoKB_cancer_gene"),
    input_table_path = req_non_existing_ann_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -OncoKB_oncogene_TSG, -EFO, -OncoKB_cancer_gene,
    EFO, OncoKB_cancer_gene))


readr::write_tsv(
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -MONDO, -OncoKB_cancer_gene, -Protein_RefSeq_ID, -RMTL),
  req_non_existing_ann_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("MONDO", "OncoKB_cancer_gene", "Protein_RefSeq_ID"),
    input_table_path = req_non_existing_ann_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -RMTL, -MONDO, -OncoKB_cancer_gene, -Protein_RefSeq_ID,
    MONDO, OncoKB_cancer_gene, Protein_RefSeq_ID))


readr::write_tsv(
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -EFO, -OncoKB_oncogene_TSG, -Gene_full_name, -RMTL),
  req_non_existing_ann_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "OncoKB_oncogene_TSG", "Gene_full_name"),
    input_table_path = req_non_existing_ann_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -RMTL, -EFO, -OncoKB_oncogene_TSG, -Gene_full_name,
    EFO, OncoKB_oncogene_TSG, Gene_full_name))


readr::write_tsv(
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -EFO, -OncoKB_oncogene_TSG, -Gene_full_name, -RMTL, -Gene_type),
  req_non_existing_ann_tbl_path)

testthat::expect_equal(
  run_cli_get_tibble(
    columns_to_add = c("EFO", "OncoKB_oncogene_TSG", "Gene_full_name",
                       "Gene_type"),
    input_table_path = req_non_existing_ann_tbl_path,
    output_table_path = annotator_cli_output_path,
    remove_input_table = TRUE),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -RMTL, -EFO, -OncoKB_oncogene_TSG, -Gene_full_name, -Gene_type,
    EFO, OncoKB_oncogene_TSG, Gene_full_name, Gene_type))

# Error on non-existing input file
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("EFO"),
    input_table_path = "non_existing_FILE.tsv",
    output_table_path = annotator_cli_output_path))

# Error on non-existing output dir
testthat::expect_warning(
  run_cli_get_tibble(
    columns_to_add = c("EFO"),
    input_table_path = working_input_tsv_path,
    output_table_path = "non_existing_DIR/test.tsv"))
