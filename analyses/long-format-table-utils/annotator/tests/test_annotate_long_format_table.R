# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_annotate_long_format_table.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
annotate_long_format_table <- import_function(
  "../annotator-api.R", "annotate_long_format_table")

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
  "test_data/test_long_format_table.tsv",
  col_types = readr::cols(.default = readr::col_character()))[]

inspected_annotated_long_format_tibble <- readr::read_tsv(
  "test_data/inspected_annotated_test_long_format_table.tsv",
  col_types = readr::cols(.default = readr::col_character()),
  na = c("NA"), quoted_na = FALSE, trim_ws = FALSE)[]

# nac means non all character. These tables can be used to test whether column
# types are changed.
long_format_tibble_nac <- readr::read_tsv(
  "test_data/test_long_format_table.tsv",
  col_types = readr::cols(
    .default = readr::col_character(),
    tpm_mean = readr::col_double()))[]

inspected_annotated_long_format_tibble_nac <- readr::read_tsv(
  "test_data/inspected_annotated_test_long_format_table.tsv",
  col_types = readr::cols(
    .default = readr::col_character(),
    tpm_mean = readr::col_double()),
  na = c("NA"), quoted_na = FALSE, trim_ws = FALSE)[]

testthat::expect_equal(
  sum(is.na(long_format_tibble_nac)), 0)

testthat::expect_equal(
  sum(is.na(inspected_annotated_long_format_tibble_nac)), 0)

testthat::expect_equal(
  sum(is.na(inspected_annotated_long_format_tibble)), 0)

testthat::expect_equal(
  sum(is.na(long_format_tibble)), 0)


# Test cases
#
# Add package prefix for auto completion purpose only
#
# Test standard usecase
testthat::expect_equal(
  annotate_long_format_table(long_format_tibble),
  inspected_annotated_long_format_tibble)

# Test replacing NAs with empty strings for **ALL COLUMNS THAT HAVE NA** in the
# output table
#
# tpm_mean has no NA, so its type should not be changed
testthat::expect_equal(
  annotate_long_format_table(long_format_tibble_nac),
  inspected_annotated_long_format_tibble_nac)

testthat::expect_false(
  is.character(long_format_tibble_nac$tpm_mean))

testthat::expect_false(
  is.character(annotate_long_format_table(long_format_tibble_nac)$tpm_mean))

testthat::expect_false(
  is.character(annotate_long_format_table(
    long_format_tibble_nac,
    replace_na_with_empty_string = FALSE)$tpm_mean))

long_format_tibble_nac_tmna <- long_format_tibble_nac
long_format_tibble_nac_tmna[2, "tpm_mean"] <- NA
testthat::expect_true(
  is.character(
    annotate_long_format_table(long_format_tibble_nac_tmna)$tpm_mean))

testthat::expect_false(
  is.character(
    annotate_long_format_table(
      long_format_tibble_nac_tmna,
      replace_na_with_empty_string = FALSE)$tpm_mean))

testthat::expect_gt(
  sum(is.na(
    annotate_long_format_table(
      long_format_tibble_nac_tmna,
      replace_na_with_empty_string = FALSE))),
  0)

testthat::expect_equal(
  sum(is.na(
    annotate_long_format_table(
      long_format_tibble_nac_tmna,
      replace_na_with_empty_string = TRUE))),
  0)

# Test annotation order
testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("RMTL", "Gene_type", "OncoKB_cancer_gene",
                       "OncoKB_oncogene_TSG", "Gene_full_name",
                       "Protein_RefSeq_ID", "EFO", "MONDO")),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean", "RMTL",
      "Gene_type", "OncoKB_cancer_gene", "OncoKB_oncogene_TSG",
      "Gene_full_name", "Protein_RefSeq_ID", "EFO", "MONDO")])

testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble, columns_to_add = c("MONDO", "RMTL", "EFO")),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "MONDO", "RMTL", "EFO")])

testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble, columns_to_add = c("RMTL", "EFO", "MONDO")),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "RMTL", "EFO", "MONDO")])

testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("RMTL", "Protein_RefSeq_ID", "Gene_full_name")),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "RMTL", "Protein_RefSeq_ID", "Gene_full_name")])

testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("OncoKB_oncogene_TSG", "OncoKB_cancer_gene", "MONDO")),
  inspected_annotated_long_format_tibble[,
    c("Gene_symbol", "Gene_Ensembl_ID", "Disease", "cohort", "tpm_mean",
      "OncoKB_oncogene_TSG", "OncoKB_cancer_gene", "MONDO")])

# Return same table if no annotation to add
testthat::expect_equal(
  annotate_long_format_table(
    long_format_tibble, columns_to_add = character(0)),
  long_format_tibble)

# Error on duplicated annotation columns
testthat::expect_error(annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("OncoKB_oncogene_TSG", "OncoKB_oncogene_TSG", "MONDO")))

testthat::expect_error(annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("OncoKB_oncogene_TSG", "MONDO", "MONDO")))

# Error on non-available annotation columns
testthat::expect_error(annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("NOT_AVAILABLE", "MONDO")))
testthat::expect_error(annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("NOT_AVAILABLE")))

# Error on invalid input table
testthat::expect_error(
  annotate_long_format_table(c(1, 2)))

testthat::expect_error(
  annotate_long_format_table(data.frame(a = c(1, 2))))

testthat::expect_error(
  annotate_long_format_table(data.table::data.table(a = c(1, 2))))

null_colname_long_format_tibble <- long_format_tibble
colnames(null_colname_long_format_tibble) <- NULL
testthat::expect_error(
  annotate_long_format_table(null_colname_long_format_tibble))

# Error on missing required columns
testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(long_format_tibble, -Gene_symbol)))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(long_format_tibble, -Gene_Ensembl_ID)))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(long_format_tibble, -Disease)))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(long_format_tibble, -Gene_symbol, -Gene_Ensembl_ID)))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(
      long_format_tibble, -Gene_symbol, -Gene_Ensembl_ID, -Disease)))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(long_format_tibble, -Gene_Ensembl_ID, -Disease)))


# Error on requiring existing annotation columns
testthat::expect_error(
  annotate_long_format_table(
    inspected_annotated_long_format_tibble))

testthat::expect_error(
  annotate_long_format_table(
    dplyr::select(inspected_annotated_long_format_tibble, -EFO),
    columns_to_add = c("EFO", "MONDO")))

# Error on duplicated annotation columns
testthat::expect_error(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c("EFO", "EFO")))

# Error on NULL annotation columns
testthat::expect_error(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = NULL))

# Error on non character annotation columns
testthat::expect_error(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c(1)))

testthat::expect_error(
  annotate_long_format_table(
    long_format_tibble,
    columns_to_add = c(TRUE)))

# No error on requiring non-existing annotation columns
#
# Relocate to last, so the order is expected. Adapted from
# https://stackoverflow.com/a/43902237/4638182. The dplyr::relocate is not
# available in the Docker image
#
# The behavior of testthat::expect_equal changed at some point. The Docker
# image/container version does not check column order, whereas the latest
# version checks.
testthat::expect_equal(
  annotate_long_format_table(
    dplyr::select(
      inspected_annotated_long_format_tibble, -EFO, -OncoKB_cancer_gene),
    columns_to_add = c("EFO", "OncoKB_cancer_gene")),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -EFO, -OncoKB_cancer_gene,
    EFO, OncoKB_cancer_gene))

testthat::expect_equal(
  annotate_long_format_table(
    dplyr::select(
      inspected_annotated_long_format_tibble,
      -EFO, -OncoKB_cancer_gene, -OncoKB_oncogene_TSG),
    columns_to_add = c("EFO", "OncoKB_cancer_gene")),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -OncoKB_oncogene_TSG, -EFO, -OncoKB_cancer_gene,
    EFO, OncoKB_cancer_gene))

testthat::expect_equal(
  annotate_long_format_table(
    dplyr::select(
      inspected_annotated_long_format_tibble,
      -MONDO, -OncoKB_cancer_gene, -Protein_RefSeq_ID, -RMTL),
    columns_to_add = c("MONDO", "OncoKB_cancer_gene", "Protein_RefSeq_ID")),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -RMTL, -MONDO, -OncoKB_cancer_gene, -Protein_RefSeq_ID,
    MONDO, OncoKB_cancer_gene, Protein_RefSeq_ID))

testthat::expect_equal(
  annotate_long_format_table(
    dplyr::select(
      inspected_annotated_long_format_tibble,
      -EFO, -OncoKB_oncogene_TSG, -Gene_full_name, -RMTL),
    columns_to_add = c("EFO", "OncoKB_oncogene_TSG", "Gene_full_name")),
  dplyr::select(
    inspected_annotated_long_format_tibble,
    -RMTL, -EFO, -OncoKB_oncogene_TSG, -Gene_full_name,
    EFO, OncoKB_oncogene_TSG, Gene_full_name))
