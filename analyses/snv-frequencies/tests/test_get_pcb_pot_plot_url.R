# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_get_pcb_pot_plot_url.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
get_pcb_pot_plot_url <- import_function(
  "../01-snv-frequencies.R", "get_pcb_pot_plot_url")

# Add package prefix for auto completion purpose only
#
# Test cases

# Return one URL for each gene, as get_pcb_pot_plot_url is used for
# dplyr::mutate
testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("BRAF", "CTNNB1"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/oncoprint?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=", c("BRAF", "CTNNB1")))

testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/oncoprint?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=BRAF"))

testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("CTNNB1"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/oncoprint?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=CTNNB1"))

# - test mutations plot
testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("BRAF", "CTNNB1"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "mutations"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/mutations?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=", c("BRAF", "CTNNB1")))

testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "mutations"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/mutations?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=BRAF"))

testthat::expect_equal(
  get_pcb_pot_plot_url(
    c("CTNNB1"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "mutations"),
  paste0(
    "https://pedcbioportal.kidsfirstdrc.org/results/mutations?",
    "cancer_study_list=ped_opentargets_2021&case_set_id=",
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "&Action=Submit&gene_list=CTNNB1"))

# Return character(0) if no gene
testthat::expect_equal(
  get_pcb_pot_plot_url(
    character(0),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"),
  character(0))
# Error on multiple plot types
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("CTNNB1"),
    "ped_opentargets_2021_pbta_medulloblastoma",
    c("mutations", "oncoprint")))
# Error on empty case_set_id
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("CTNNB1", "CTNNB1"), character(0),
    "mutations"))

testthat::expect_error(
  get_pcb_pot_plot_url(
    c("CTNNB1", "CTNNB1"), character(0),
    "oncoprint"))
# Error on multiple case_set_ids
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("CTNNB1", "CTNNB1"),
    c("ped_opentargets_2021_pbta_medulloblastoma",
      "ped_opentargets_2021_pbta_medulloblastoma"),
    "oncoprint"))

testthat::expect_error(
  get_pcb_pot_plot_url(
    c("CTNNB1", "CTNNB1"),
    c("ped_opentargets_2021_pbta_medulloblastoma",
      "ped_opentargets_2021_pbta_medulloblastoma"),
    "mutations"))
# Error on NULL
testthat::expect_error(
  get_pcb_pot_plot_url(
    NULL,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    NULL,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NULL))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NULL,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NULL))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NULL,
    NULL,
    NULL))

# Error on NA
testthat::expect_error(
  get_pcb_pot_plot_url(
    NA_character_,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    NA_character_,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NA_character_))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NA_character_,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NA_character_))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NA_character_,
    NA_character_,
    NA_character_))

testthat::expect_error(
  get_pcb_pot_plot_url(
    NaN,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    NaN,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NaN))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NaN,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    NaN))
testthat::expect_error(
  get_pcb_pot_plot_url(
    NaN,
    NaN,
    NaN))
# Error on non-character
testthat::expect_error(
  get_pcb_pot_plot_url(
    1,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    1,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    1))
testthat::expect_error(
  get_pcb_pot_plot_url(
    1,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    1))
testthat::expect_error(
  get_pcb_pot_plot_url(
    1,
    1,
    1))

testthat::expect_error(
  get_pcb_pot_plot_url(
    TRUE,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    TRUE,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    TRUE))
testthat::expect_error(
  get_pcb_pot_plot_url(
    TRUE,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    TRUE))
testthat::expect_error(
  get_pcb_pot_plot_url(
    TRUE,
    TRUE,
    TRUE))

testthat::expect_error(
  get_pcb_pot_plot_url(
    Inf,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    Inf,
    "oncoprint"))
testthat::expect_error(
  get_pcb_pot_plot_url(
    c("BRAF"),
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    Inf))
testthat::expect_error(
  get_pcb_pot_plot_url(
    Inf,
    "ped_opentargets_2021_high-grade_glioma_astrocytoma",
    Inf))
testthat::expect_error(
  get_pcb_pot_plot_url(
    Inf,
    Inf,
    Inf))
