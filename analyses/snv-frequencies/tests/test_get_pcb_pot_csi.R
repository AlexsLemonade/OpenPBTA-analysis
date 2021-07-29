# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_get_pcb_pot_csi.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
get_pcb_pot_csi <- import_function(
  "../01-snv-frequencies.R", "get_pcb_pot_csi")

# Add package prefix for auto completion purpose only
#
# Test cases
testthat::expect_equal(
  get_pcb_pot_csi(c("High-grade glioma/astrocytoma"), c("PBTA", "GMKF")),
  "ped_opentargets_2021_high-grade_glioma_astrocytoma")
testthat::expect_equal(
  get_pcb_pot_csi(c("High-grade glioma/astrocytoma"), c("GMKF")),
  "ped_opentargets_2021_gmkf_high-grade_glioma_astrocytoma")
testthat::expect_equal(
  get_pcb_pot_csi(c("High-grade glioma/astrocytoma"), c("NULL")),
  "ped_opentargets_2021_null_high-grade_glioma_astrocytoma")
testthat::expect_equal(
  get_pcb_pot_csi(c("High-grade glioma/;=!@#astrocytoma"), c("GMKF")),
  "ped_opentargets_2021_gmkf_high-grade_glioma______astrocytoma")
testthat::expect_equal(
  get_pcb_pot_csi(c("High-grade glioma123456astrocytoma"), c("GMKF")),
  "ped_opentargets_2021_gmkf_high-grade_glioma123456astrocytoma")
testthat::expect_equal(
  get_pcb_pot_csi(c("Medulloblastoma"), c("PBTA")),
  "ped_opentargets_2021_pbta_medulloblastoma")

# Error on empty characters
testthat::expect_error(get_pcb_pot_csi(character(0), character(0)))
testthat::expect_error(get_pcb_pot_csi(c(), c()))
testthat::expect_error(get_pcb_pot_csi(character(0), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", character(0)))
testthat::expect_error(get_pcb_pot_csi(c(), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c()))
# Error on non-character input
testthat::expect_error(get_pcb_pot_csi(c(1), c(1)))
testthat::expect_error(get_pcb_pot_csi(c(1, 2), c(1, 2)))
testthat::expect_error(get_pcb_pot_csi(c(1), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(1)))
testthat::expect_error(get_pcb_pot_csi(c(TRUE), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(TRUE)))
testthat::expect_error(get_pcb_pot_csi(c(Inf), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(Inf)))
testthat::expect_error(get_pcb_pot_csi(c(-Inf), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(-Inf)))
# Error on NA
testthat::expect_error(get_pcb_pot_csi(c(NA_character_), c(NA_character_)))
testthat::expect_error(get_pcb_pot_csi(c(NaN), c(NaN)))
testthat::expect_error(get_pcb_pot_csi(c(NA_character_), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(NA_character_)))
testthat::expect_error(get_pcb_pot_csi(c(NA), "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", c(NA)))
# Error on NULL
testthat::expect_error(get_pcb_pot_csi(NULL, "PBTA"))
testthat::expect_error(get_pcb_pot_csi("Medulloblastoma", NULL))
testthat::expect_error(get_pcb_pot_csi(NULL, NULL))
# Error on more than one cancer_groups
testthat::expect_error(
  get_pcb_pot_csi(c("Medulloblastoma", "High-grade glioma/astrocytoma"),
                  c("PBTA")))
testthat::expect_error(
  get_pcb_pot_csi(c("Medulloblastoma", "Medulloblastoma"),
                  c("PBTA")))
