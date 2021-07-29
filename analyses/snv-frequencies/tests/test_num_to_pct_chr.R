# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_num_to_pct_chr.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
num_to_pct_chr <- import_function(
  "../01-snv-frequencies.R", "num_to_pct_chr")

# Add package prefix for auto completion purpose only
#
# Test cases
testthat::expect_equal(
  num_to_pct_chr(numeric(0)), character(0))
testthat::expect_equal(
  num_to_pct_chr(c(1, -1.1, 0, 0.1, 0/0)),
  c("100.00%", "-110.00%", "0.00%", "10.00%", ""))
testthat::expect_equal(
  num_to_pct_chr(c(1, -1.1, 0, 0.1, NaN)),
  c("100.00%", "-110.00%", "0.00%", "10.00%", ""))

# Error on non-numeric empty vector
testthat::expect_error(num_to_pct_chr(character(0)))
testthat::expect_error(num_to_pct_chr(logical(0)))
# Error on NULL
testthat::expect_error(num_to_pct_chr(NULL))
# Error on Inf
testthat::expect_error(num_to_pct_chr(c(1, -1.1, 0, 0.1, 1/0)))
testthat::expect_error(num_to_pct_chr(c(1, Inf)))
# Error on -Inf
testthat::expect_error(num_to_pct_chr(c(1, -1.1, 0, 0.1, -1/0)))
testthat::expect_error(num_to_pct_chr(c(1, -Inf)))
# Error on NA
testthat::expect_error(num_to_pct_chr(c(1, NA)))
# Error on NA and Inf
testthat::expect_error(num_to_pct_chr(c(1, -1.1, 0, 0.1, NA, 0.1/0)))
