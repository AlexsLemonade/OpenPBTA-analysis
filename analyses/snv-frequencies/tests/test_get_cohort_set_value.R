# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_get_cohort_set_value.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
get_cohort_set_value <- import_function(
  "../01-snv-frequencies.R", "get_cohort_set_value")

# Add package prefix for auto completion purpose only
#
# Test cases
# Return a single cohort if there is only one cohort
testthat::expect_equal(get_cohort_set_value(c("a")), "a")
testthat::expect_equal(get_cohort_set_value(c("ABC")), "ABC")
testthat::expect_equal(get_cohort_set_value(c("AbC")), "AbC")
testthat::expect_equal(get_cohort_set_value(c("1")), "1")
# - retain special chars
testthat::expect_equal(
  get_cohort_set_value(c("1/2/3;'!@#$%^&*(()-<-<>")),
  "1/2/3;'!@#$%^&*(()-<-<>")
# Return "all_cohorts" if there are more than one cohorts
testthat::expect_equal(get_cohort_set_value(c("a", "b")), "all_cohorts")
testthat::expect_equal(get_cohort_set_value(c("A", "B", "C")), "all_cohorts")
testthat::expect_equal(get_cohort_set_value(c("a", "b", "c")), "all_cohorts")
# - handle special chars
testthat::expect_equal(
  get_cohort_set_value(c("1/2/3;'!@#$%^&*(()-", "*&^%&#*(!@\"';:>")),
  "all_cohorts")
# - case sensitive
testthat::expect_equal(get_cohort_set_value(c("a", "A")), "all_cohorts")
testthat::expect_equal(get_cohort_set_value(c("abc", "abC")), "all_cohorts")

# Error on duplicated cohort values
testthat::expect_error(get_cohort_set_value("a", "a"))
testthat::expect_error(get_cohort_set_value("A", "A"))
# Error on NA value
testthat::expect_error(get_cohort_set_value(c("a", NA)))
testthat::expect_error(get_cohort_set_value(c(NA_character_)))
testthat::expect_error(get_cohort_set_value(c(NA_integer_)))
testthat::expect_error(get_cohort_set_value(NA))
# Error on NULL
testthat::expect_error(get_cohort_set_value(NULL))
# Error on empty character cohort vector
testthat::expect_error(get_cohort_set_value(character(0)))
# Error on empty numeric cohort vector
testthat::expect_error(get_cohort_set_value(numeric(0)))
# Error on non-empty non-character cohort vector
testthat::expect_error(get_cohort_set_value(c(1)))
testthat::expect_error(get_cohort_set_value(c(TRUE)))
testthat::expect_error(get_cohort_set_value(c(1, 2)))
