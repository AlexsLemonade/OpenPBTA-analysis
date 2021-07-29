# The working directory is the directory that contains this test R file, if this
# file is executed by test_dir
#
# testthat package is loaded, if this file is executed by test_dir
context("tests/test_get_opr_mut_freq_tbl.R")
# import_function is defined in tests/helper_import_function.R and tested in
# annotator/tests/test_helper_import_function.R
#
# testthat::test_dir("tests") runs all `helper*R` files under the `tests`
# directory before running the `test*R` files
get_opr_mut_freq_tbl <- import_function(
  "../01-snv-frequencies.R", "get_opr_mut_freq_tbl")

# load dependency functions
library(dplyr)
num_to_pct_chr <- import_function(
  "../01-snv-frequencies.R", "num_to_pct_chr")

# Craft some data for testing
test_maf_tbl <- tibble::tribble(
  ~Kids_First_Biospecimen_ID, ~Gene, ~Var,
  "BS_1", "GeneA", "Var1",
  "BS_2", "GeneA", "Var1",
  "BS_2", "GeneA", "Var2",
  "BS_2", "GeneB", "Var3",
  "BS_3", "GeneB", "Var3",
  "BS_3", "GeneC", "Var4",
  "BS_3", "GeneC", "Var5",
  "BS_4", "GeneB", "Var3",
  "BS_5", "GeneC", "Var5",
)

test_overall_histology_tbl <- tibble::tribble(
  ~Kids_First_Participant_ID, ~Kids_First_Biospecimen_ID,
  "PT_1", "BS_1",
  "PT_1", "BS_2",
  "PT_2", "BS_3",
  "PT_2", "BS_4",
  "PT_3", "BS_5",
)

# indp means independent sample
#
# computed using each cohort or all cohorts procedures
test_primary_indp_each_tbl <- tibble::tribble(
  ~Kids_First_Participant_ID, ~Kids_First_Biospecimen_ID,
  "PT_1", "BS_1",
  "PT_1", "BS_2",
  "PT_2", "BS_3",
  "PT_3", "BS_5",
)

test_primary_indp_all_tbl <- tibble::tribble(
  ~Kids_First_Participant_ID, ~Kids_First_Biospecimen_ID,
  "PT_1", "BS_1",
  "PT_2", "BS_3",
  "PT_3", "BS_5",
)

test_relapse_indp_all_tbl <- tibble::tribble(
  ~Kids_First_Participant_ID, ~Kids_First_Biospecimen_ID,
  "PT_2", "BS_4",
  "PT_3", "BS_5",
)

test_relapse_indp_each_tbl <- tibble::tribble(
  ~Kids_First_Participant_ID, ~Kids_First_Biospecimen_ID,
  "PT_2", "BS_4",
  "PT_3", "BS_5",
)

test_each_indp_gene_mut_freq_tbl <- tibble::tibble(
  Gene = c("GeneA", "GeneB", "GeneC"),
  # number of patients with the mutation
  Total_mutations = as.integer(c(1, 2, 2)),
  Patients_in_dataset = as.integer(3),
  Total_mutations_Over_Patients_in_dataset = c("1/3", "2/3", "2/3"),
  Frequency_in_overall_dataset = num_to_pct_chr(c(1/3, 2/3, 2/3)),
  # number of primary samples with the mutation
  Total_primary_tumors_mutated = as.integer(c(2, 2, 2)),
  Total_relapse_tumors_mutated = as.integer(c(0, 1, 1)),
  Primary_tumors_in_dataset = as.integer(4),
  Relapse_tumors_in_dataset = as.integer(2),
  Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = c(
    "2/4", "2/4", "2/4"
  ),
  Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = c(
    "0/2", "1/2", "1/2"
  ),
  Frequency_in_primary_tumors = num_to_pct_chr(c(2/4, 2/4, 2/4)),
  Frequency_in_relapse_tumors = num_to_pct_chr(c(0/2, 1/2, 1/2))
)

test_all_indp_gene_mut_freq_tbl <- tibble::tibble(
  Gene = c("GeneA", "GeneB", "GeneC"),
  # number of patients with the mutation
  Total_mutations = as.integer(c(1, 2, 2)),
  Patients_in_dataset = as.integer(3),
  Total_mutations_Over_Patients_in_dataset = c("1/3", "2/3", "2/3"),
  Frequency_in_overall_dataset = num_to_pct_chr(c(1/3, 2/3, 2/3)),
  # number of primary samples with the mutation
  Total_primary_tumors_mutated = as.integer(c(1, 1, 2)),
  Total_relapse_tumors_mutated = as.integer(c(0, 1, 1)),
  Primary_tumors_in_dataset = as.integer(3),
  Relapse_tumors_in_dataset = as.integer(2),
  Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = c(
    "1/3", "1/3", "2/3"
  ),
  Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = c(
    "0/2", "1/2", "1/2"
  ),
  Frequency_in_primary_tumors = num_to_pct_chr(c(1/3, 1/3, 2/3)),
  Frequency_in_relapse_tumors = num_to_pct_chr(c(0/2, 1/2, 1/2))
)

test_each_indp_var_mut_freq_tbl <- tibble::tibble(
  Var = c("Var1", "Var2", "Var3", "Var4", "Var5"),
  # number of patients with the mutation
  Total_mutations = as.integer(c(1, 1, 2, 1, 2)),
  Patients_in_dataset = as.integer(3),
  Total_mutations_Over_Patients_in_dataset = c(
    "1/3", "1/3", "2/3", "1/3", "2/3"
  ),
  Frequency_in_overall_dataset = num_to_pct_chr(c(1/3, 1/3, 2/3, 1/3, 2/3)),
  # number of primary samples with the mutation
  Total_primary_tumors_mutated = as.integer(c(2, 1, 2, 1, 2)),
  Total_relapse_tumors_mutated = as.integer(c(0, 0, 1, 0, 1)),
  Primary_tumors_in_dataset = as.integer(4),
  Relapse_tumors_in_dataset = as.integer(2),
  Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = c(
    "2/4", "1/4", "2/4", "1/4", "2/4"
  ),
  Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = c(
    "0/2", "0/2", "1/2", "0/2", "1/2"
  ),
  Frequency_in_primary_tumors = num_to_pct_chr(c(2/4, 1/4, 2/4, 1/4, 2/4)),
  Frequency_in_relapse_tumors = num_to_pct_chr(c(0/2, 0/2, 1/2, 0/2, 1/2))
)

test_all_indp_var_mut_freq_tbl <- tibble::tibble(
  Var = c("Var1", "Var2", "Var3", "Var4", "Var5"),
  # number of patients with the mutation
  Total_mutations = as.integer(c(1, 1, 2, 1, 2)),
  Patients_in_dataset = as.integer(3),
  Total_mutations_Over_Patients_in_dataset = c(
    "1/3", "1/3", "2/3", "1/3", "2/3"
  ),
  Frequency_in_overall_dataset = num_to_pct_chr(c(1/3, 1/3, 2/3, 1/3, 2/3)),
  # number of primary samples with the mutation
  Total_primary_tumors_mutated = as.integer(c(1, 0, 1, 1, 2)),
  Total_relapse_tumors_mutated = as.integer(c(0, 0, 1, 0, 1)),
  Primary_tumors_in_dataset = as.integer(3),
  Relapse_tumors_in_dataset = as.integer(2),
  Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = c(
    "1/3", "0/3", "1/3", "1/3", "2/3"
  ),
  Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = c(
    "0/2", "0/2", "1/2", "0/2", "1/2"
  ),
  Frequency_in_primary_tumors = num_to_pct_chr(c(1/3, 0/3, 1/3, 1/3, 2/3)),
  Frequency_in_relapse_tumors = num_to_pct_chr(c(0/2, 0/2, 1/2, 0/2, 1/2))
)

# Add package prefix for auto completion purpose only
#
# Test cases
testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl),
  test_all_indp_gene_mut_freq_tbl,
)

testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    test_primary_indp_each_tbl, test_relapse_indp_each_tbl),
  test_each_indp_gene_mut_freq_tbl
)

testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Var", test_overall_histology_tbl,
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl),
  test_all_indp_var_mut_freq_tbl
)

testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Var", test_overall_histology_tbl,
    test_primary_indp_each_tbl, test_relapse_indp_each_tbl),
  test_each_indp_var_mut_freq_tbl
)
# Test 0/0
testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    test_relapse_indp_all_tbl),
  test_all_indp_gene_mut_freq_tbl %>% mutate(
    Total_primary_tumors_mutated = as.integer(0),
    Primary_tumors_in_dataset = as.integer(0),
    Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = "0/0",
    Frequency_in_primary_tumors = num_to_pct_chr(0/0)
  ),
)

testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    test_relapse_indp_each_tbl),
  test_each_indp_gene_mut_freq_tbl %>% mutate(
    Total_primary_tumors_mutated = as.integer(0),
    Primary_tumors_in_dataset = as.integer(0),
    Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = "0/0",
    Frequency_in_primary_tumors = num_to_pct_chr(0/0)
  ),
)


testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0))),
  test_all_indp_gene_mut_freq_tbl %>% mutate(
    Total_relapse_tumors_mutated = as.integer(0),
    Relapse_tumors_in_dataset = as.integer(0),
    Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = "0/0",
    Frequency_in_relapse_tumors = num_to_pct_chr(0/0)
  ),
)



testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    test_primary_indp_each_tbl,
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0))),
  test_each_indp_gene_mut_freq_tbl %>% mutate(
    Total_relapse_tumors_mutated = as.integer(0),
    Relapse_tumors_in_dataset = as.integer(0),
    Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = "0/0",
    Frequency_in_relapse_tumors = num_to_pct_chr(0/0)
  ),
)


testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene", test_overall_histology_tbl,
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0))),
  test_each_indp_gene_mut_freq_tbl %>% mutate(
    Total_primary_tumors_mutated = as.integer(0),
    Primary_tumors_in_dataset = as.integer(0),
    Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = "0/0",
    Frequency_in_primary_tumors = num_to_pct_chr(0/0),
    Total_relapse_tumors_mutated = as.integer(0),
    Relapse_tumors_in_dataset = as.integer(0),
    Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = "0/0",
    Frequency_in_relapse_tumors = num_to_pct_chr(0/0)
  ),
)

# Test all empty case
testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl[NULL, ], "Gene",
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0))),
  test_each_indp_gene_mut_freq_tbl[NULL, ],
)

testthat::expect_equal(
  get_opr_mut_freq_tbl(
    test_maf_tbl[NULL, ], "Var",
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0))),
  test_each_indp_var_mut_freq_tbl[NULL, ],
)

# Error if test_maf_tbl Kids_First_Biospecimen_ID not in overall_histology_df
testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Gene",
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)))
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "Var",
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)),
    tibble::tibble(
      Kids_First_Biospecimen_ID = character(0),
      Kids_First_Participant_ID = character(0)))
)

# Error on missing group column
testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl, "UNAVAILABLE", test_overall_histology_tbl,
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl)
)

# Error on missing Kids_First_Biospecimen_ID column
testthat::expect_error(
  get_opr_mut_freq_tbl(
    select(test_maf_tbl, -Kids_First_Biospecimen_ID),
    "Gene", test_overall_histology_tbl,
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    select(test_maf_tbl, -Kids_First_Biospecimen_ID),
    "Var", test_overall_histology_tbl,
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", select(test_overall_histology_tbl, -Kids_First_Biospecimen_ID),
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", test_overall_histology_tbl,
    select(test_primary_indp_all_tbl, -Kids_First_Biospecimen_ID),
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    select(test_relapse_indp_all_tbl, -Kids_First_Biospecimen_ID))
)


# Error on missing Kids_First_Participant_ID column
testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", select(test_overall_histology_tbl, -Kids_First_Participant_ID),
    test_primary_indp_all_tbl, test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", test_overall_histology_tbl,
    select(test_primary_indp_all_tbl, -Kids_First_Participant_ID),
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene", test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    select(test_relapse_indp_all_tbl, -Kids_First_Participant_ID))
)

# Error on NA in Kids_First_Biospecimen_ID or Kids_First_Participant_ID column
testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    mutate(
      test_overall_histology_tbl,
      Kids_First_Participant_ID = if_else(
        stringr::str_detect(Kids_First_Participant_ID, "2"),
        NA, Kids_First_Participant_ID)),
    test_primary_indp_all_tbl,
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    mutate(
      test_primary_indp_all_tbl,
      Kids_First_Participant_ID = if_else(
        stringr::str_detect(Kids_First_Participant_ID, "2"),
        NA, Kids_First_Participant_ID)),
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    mutate(
      test_relapse_indp_all_tbl,
      Kids_First_Participant_ID = if_else(
        stringr::str_detect(Kids_First_Participant_ID, "2"),
        NA, Kids_First_Participant_ID)))
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    mutate(
      test_overall_histology_tbl,
      Kids_First_Biospecimen_ID = if_else(
        stringr::str_detect(Kids_First_Biospecimen_ID, "2"),
        NA, Kids_First_Biospecimen_ID)),
    test_primary_indp_all_tbl,
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    mutate(
      test_primary_indp_all_tbl,
      Kids_First_Biospecimen_ID = if_else(
        stringr::str_detect(Kids_First_Biospecimen_ID, "2"),
        NA, Kids_First_Biospecimen_ID)),
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    mutate(
      test_relapse_indp_all_tbl,
      Kids_First_Biospecimen_ID = if_else(
        stringr::str_detect(Kids_First_Biospecimen_ID, "2"),
        NA, Kids_First_Biospecimen_ID)))
)


testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    mutate(
      test_overall_histology_tbl,
      Kids_First_Biospecimen_ID = NA),
    test_primary_indp_all_tbl,
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    mutate(
      test_primary_indp_all_tbl,
      Kids_First_Biospecimen_ID = NA),
    test_relapse_indp_all_tbl)
)

testthat::expect_error(
  get_opr_mut_freq_tbl(
    test_maf_tbl,
    "Gene",
    test_overall_histology_tbl,
    test_primary_indp_all_tbl,
    mutate(
      test_relapse_indp_all_tbl,
      Kids_First_Biospecimen_ID = NA))
)
