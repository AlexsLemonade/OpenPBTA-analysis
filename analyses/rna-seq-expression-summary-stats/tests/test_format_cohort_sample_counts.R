# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_format_cohort_sample_counts.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
format_cohort_sample_counts <- import_function(
  "../01-tpm-summary-stats.R", "format_cohort_sample_counts")
library(tibble)
library(dplyr)

# Add package prefix for auto completion purpose only
#
# Test cases
testthat::expect_equal(
  format_cohort_sample_counts(c("b", "a", "a", "c", "c")),
  "a_2_samples&b_1_samples&c_2_samples")
testthat::expect_equal(
  format_cohort_sample_counts(c("a", "a")),
  "a_2_samples")
testthat::expect_equal(
  format_cohort_sample_counts(c("a")),
  "a_1_samples")
# Return empty string if input is empty
testthat::expect_equal(
  format_cohort_sample_counts(character(0)),
  "")
# Case sensitive
testthat::expect_equal(
  format_cohort_sample_counts(c("b", "A", "a", "a", "B", "c", "c")),
  "a_2_samples&A_1_samples&b_1_samples&B_1_samples&c_2_samples")
# Special chars
testthat::expect_equal(
  format_cohort_sample_counts(
    c("b;", "b", "B", "A/", "A/", "a^", "a^", "a^", "a%", "B!", "c-", "c_")),
    paste0(
      "A/_2_samples&a%_1_samples&a^_3_samples",
      "&b_1_samples&B_1_samples&b;_1_samples&B!_1_samples",
      "&c__1_samples&c-_1_samples"))

# Error on NA
testthat::expect_error(format_cohort_sample_counts(c("a", NA)))
testthat::expect_error(format_cohort_sample_counts(c(NA_character_)))
testthat::expect_error(
  format_cohort_sample_counts(c(NA_character_, NA_character_)))
testthat::expect_error(format_cohort_sample_counts(c(NaN)))
# Error on non-character
testthat::expect_error(format_cohort_sample_counts(NULL))
testthat::expect_error(format_cohort_sample_counts(NA_integer_))
testthat::expect_error(format_cohort_sample_counts(1))
testthat::expect_error(format_cohort_sample_counts(c(1,2)))
testthat::expect_error(format_cohort_sample_counts(c(TRUE, FALSE)))
