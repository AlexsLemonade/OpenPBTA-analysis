# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_helper_import_function.R")
# import_function is defined in tests/helper_import_function.R and tested here
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files

# Test cases
#
# Standard use case: import a function that is defined only once and is not
# nested in other expressions
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

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_assign_defined_ret2")(),
  2)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret2")(),
  2)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret3")(),
  3)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret4")(),
  4)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret5")(),
  5)

expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret5")(),
  5)

# single_lhs_ld_rhs_defined_ret6 is defined to call
# single_lhs_ld_rhs_defined_ret5
#
# which is not defined here yet
expect_error(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret6")())
# Import single_lhs_ld_rhs_defined_ret5
single_lhs_ld_rhs_defined_ret5 <- import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret5")
expect_equal(
  import_function(
    "test_data/test_import_function_non_empty.R",
    "single_lhs_ld_rhs_defined_ret6")(),
  6)

# Error on importing from empty file
expect_error(import_function("test_data/test_import_function_empty.R", "foo"))

# Error on importing functions defined multiple times
expect_error(import_function(
  "test_data/test_import_function_non_empty.R", "multi_defined"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R", "multi_defined_multi_ways1"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R", "multi_defined_multi_ways2"))

# Error on importing functions that is nested in other expressions
expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_single_lhs_ld_rhs_defined_ret1"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_single_lhs_ld_rhs_defined_ret2"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_single_lhs_ld_rhs_defined_ret3"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_assign_defined_ret3"))


expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_single_ld_lhs_rhs_defined_ret4"))

expect_error(import_function(
  "test_data/test_import_function_non_empty.R",
  "nested_single_ld_lhs_rhs_defined_ret5"))
