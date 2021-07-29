# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_collapse_rp_lists.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
collapse_rp_lists <- import_function(
  "../download-annotation-data.R", "collapse_rp_lists")

# Test cases:
# Return NA on all NULL values
expect_equal(collapse_rp_lists(list(NULL, NULL, NULL)), NA_character_)
# Return NA on empty list
expect_equal(collapse_rp_lists(list()), NA_character_)
# Remove duplicates
expect_equal(
  collapse_rp_lists(list(c("NP_1", "NP_2", "NP_1", "NP_3"))),
  "NP_1,NP_2,NP_3")
# Only keep IDs with NP_ prefix
expect_equal(
  collapse_rp_lists(list(c("NX_5", "NP3", "NP", "NA"), c("NP_1", "NP_3"))),
  "NP_1,NP_3")
# Only keep non NULL values
expect_equal(
  collapse_rp_lists(list(c("NP_1"), NULL, c("NP_1", "NP_2"))),
  "NP_1,NP_2")
# Return values are sorted
expect_equal(
  collapse_rp_lists(list(c("NP_3", "NP_1"))),
  "NP_1,NP_3")
# Only keep non NA values
expect_equal(
  collapse_rp_lists(list(c("NP_1", "NP_2"), c("NP_1", "NP_3", NA),
                         c(NA_character_, NA_character_), c(NA_character_),
                         character(0))),
  "NP_1,NP_2,NP_3")
