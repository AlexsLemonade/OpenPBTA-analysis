# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_collapse_rp_lists.R")

expect_error(import_function("test_data/test_import_function_empty.R", "foo"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R", "multi_defined"))

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret1")(),
  1)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_assign_defined_ret2")(),
  2)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_ld_lhs_rhs_defined_ret3")(),
  3)
