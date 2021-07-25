# The working directory is the directory that contains this test R file, if this
# file is executed by testthat::test_dir
#
# testthat package is loaded, if this file is executed by testthat::test_dir
context("tests/test_collapse_name_vec.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
collapse_name_vec <- import_function(
  "../download-annotation-data.R", "collapse_name_vec")

# Test cases

# Regular case
expect_equal(collapse_name_vec(c("a", "b")), "a,b")
# NA is removed
expect_equal(collapse_name_vec(c("a", "b", NA)), "a,b")
# Duplicates are removed
expect_equal(collapse_name_vec(c("a", "a", "a", "b", "d")), "a,b,d")
# Values are sorted
expect_equal(collapse_name_vec(c("d", "a", "a", "b", "d")), "a,b,d")
# Value is NA if only value is NA
expect_equal(collapse_name_vec(NA_character_), NA_character_)
# Value is a single NA
expect_equal(
  collapse_name_vec(c(NA_character_, NA_character_, NA_character_)),
  NA_character_)
# Value is NA if input is empty
expect_equal(collapse_name_vec(character(0)), NA_character_)
# Error on numeric input
expect_error(collapse_name_vec(c(1)))
