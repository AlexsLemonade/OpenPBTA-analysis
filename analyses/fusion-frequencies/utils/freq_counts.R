suppressPackageStartupMessages(library(tidyverse))


# Format numbers to percentage characters
# Adapted from @Richie Cotton's answer at
# https://stackoverflow.com/a/7146270/4638182
num_to_pct_chr <- function(x, digits = 2, format = "f", ...) {
  stopifnot(!is.null(x))
  if(length(x) == 0) {
    stopifnot(is.numeric(x))
    return(character(0))
  }
  
  fx <- paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  # If there is any abnormal values other than NaN, check code correctness and
  # handle them.
  stopifnot(identical(sum(is.infinite(x)), as.integer(0)))
  stopifnot(identical(fx == 'NaN%', is.na(x)))
  # replace 'NaN%' with ''
  rfx <- replace(fx, fx == 'NaN%', '')
  return(rfx)
}
# # test cases
# num_to_pct_chr(numeric(0))
# num_to_pct_chr(c(1, -1.1, 0, 0.1, 0/0))
# num_to_pct_chr(c(1, -1.1, 0, 0.1, NaN))
# # following cases should fail
# num_to_pct_chr(character(0))
# num_to_pct_chr(c())
# num_to_pct_chr(NULL)
# num_to_pct_chr(c(1, -1.1, 0, 0.1, 1/0))
# num_to_pct_chr(c(1, -1.1, 0, 0.1, -1/0))
# num_to_pct_chr(c(1, -1.1, 0, 0.1, NA))
# num_to_pct_chr(c(1, -1.1, 0, 0.1, NA, 0.1/0))



# Collapse a list of refseq.protein character vectors from
# mygene.info query results
collapse_rp_lists <- function(xl) {
  fxl <- discard(xl, is.null)
  # assert all characters
  sapply(fxl, function(x) {
    stopifnot(is.character(x))
  })
  fxv <- unique(purrr::reduce(fxl, c, .init = character(0)))
  np_fxv <- fxv[str_detect(fxv, '^NP_')]
  c_np_fxv <- paste(np_fxv, collapse = ',')
  return(c_np_fxv)
}
# # test cases
# collapse_rp_lists(list(NULL, NULL, NULL))
# collapse_rp_lists(list())
# collapse_rp_lists(list(c('NP_1', 'NP_2'), c('NP_1', 'NP_3')))
# collapse_rp_lists(list(c('NP_1', 'NP_2'), NULL, c('NP_1', 'NP_3')))
# collapse_rp_lists(list(c('NP_1', 'NP_2')))
#
# tmp <- as_tibble(
#   mg_qres_list$response[, c('query', 'notfound', 'name', 'refseq.protein')]) %>%
#   replace_na(list(notfound = FALSE)) %>%
#   filter(!notfound) %>%
#   filter(query %in% c('ENSG00000229425', 'ENSG00000171988')) %>%
#   group_by(query) %>%
#   summarise(name = paste(unique(name), collapse = ', '),
#             refseq_protein = collapse_rp_lists(refseq.protein))



# Generate a summary table of cancer_group and cohort
#
# Args:
# - histology_df: a tibble of histology information. Must contain the following
# fields: Kids_First_Biospecimen_ID, cancer_group, cohort.
#
# Returns a tibble. For each cancer_group, create a record of each cohort and
# a comma-separated-list of all cohorts. Combine and remove duplicates.
get_cg_cs_tbl <- function(histology_df) {
  fh_df <- histology_df %>%
    filter(!is.na(Kids_First_Biospecimen_ID),
           !is.na(cancer_group),
           !is.na(cohort)) %>%
    select(Kids_First_Biospecimen_ID, cancer_group, cohort) %>%
    distinct()
  
  fh_1cg_1c_df <- fh_df %>%
    group_by(cancer_group, cohort) %>%
    summarise(n_samples = n()) %>%
    mutate(cohort_list = vapply(cohort, list, list(character(1)))) %>%
    ungroup()
  
  fh_1cg_all_cs_df <- fh_df %>%
    group_by(cancer_group) %>%
    summarise(n_samples = n(),
              cohort_list = list(unique(cohort)),
              cohort = paste(unique(cohort), collapse = '&'))
  
  fh_cgc_df <- bind_rows(fh_1cg_1c_df, fh_1cg_all_cs_df) %>%
    distinct(cancer_group, cohort, n_samples, .keep_all = TRUE)
  
  # return(list(fh_1cg_1c_df, fh_1cg_all_cs_df, fh_cgc_df))
  return(fh_cgc_df)
}



# Generate mutation frequency table for (cancer_group, cohort)
#
# Args:
# - alt_df : Dataframe with Kids_First_Biospecimen_ID and a Alt_ID. For SNV  `Chromosome_Start_Position_Reference_Allele_Tumor_Seq_Allele2`, for Fusion ` FusionName_Fusion_Type` and for CNV it's `Gene_Type` where Type is deep deletion/loss/gain/amplification.
# - overall_histology_df: the histology tibble that contains all samples. Must
#   contain the following fields: Kids_First_Biospecimen_ID,
#   Kids_First_Participant_ID, cancer_group, and cohort. This
#   is used for computing the following columns: Total_alterations,
#   Patients_in_dataset, Total_alterations_Over_Patients_in_dataset and
#   Frequency_in_overall_dataset.
# - primary_histology_df: the histology tibble that contains primary tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - relapse_histology_df: the histology tibble that contains relapse tumor
#   samples.
#   Must contain the Kids_First_Biospecimen_ID field.
# - ss_cancer_group: a character value of the cancer group to compute for.
# - ss_cohorts: a vector of character values of the cohorts to compute for.
#
# Returns a tibble with additional frequency columns.
get_cg_ch_mut_freq_tbl <- function(alt_df, overall_histology_df,
                                   primary_histology_df, relapse_histology_df,
                                   ss_cancer_group, ss_cohorts) {
  stopifnot(identical(length(ss_cancer_group), as.integer(1)))
  stopifnot(identical(sum(is.na(ss_cancer_group)), as.integer(0)))
  stopifnot(!is.null(length(ss_cohorts)))
  stopifnot(identical(length(ss_cohorts), length(unique(ss_cohorts))))
  stopifnot(identical(sum(is.na(ss_cohorts)), as.integer(0)))
  stopifnot(length(ss_cohorts) >= 1)
  # ss = subset
  ss_htl_df <- overall_histology_df %>%
    filter(cancer_group == ss_cancer_group,
           cohort %in% ss_cohorts)
  
  bpid_ss_htl_df <- ss_htl_df %>%
    select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID)
  
  ss_n_patients <- length(unique(bpid_ss_htl_df$Kids_First_Participant_ID))
  
  ss_primary_kfbids <- primary_histology_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             bpid_ss_htl_df$Kids_First_Biospecimen_ID) %>%
    pull(Kids_First_Biospecimen_ID)
  
  ss_n_primary_tumors <- length(ss_primary_kfbids)
  
  ss_relapse_kfbids <- relapse_histology_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             bpid_ss_htl_df$Kids_First_Biospecimen_ID) %>%
    pull(Kids_First_Biospecimen_ID)
  
  ss_n_relapse_tumors <- length(ss_relapse_kfbids)
  
  
  ss_alt_df <- alt_df %>%
    filter(Kids_First_Biospecimen_ID %in% ss_htl_df$Kids_First_Biospecimen_ID)
  
  # gc(reset = TRUE) If one Alt_ID has more than one records in one
  # Kids_First_Biospecimen_ID, treat multiple records as a single mutation. This
  # is reasonable when computing mutation frequencies.
  #
  # Patient level frequencies: for each variant, check how many unique
  # Kids_First_Participant_IDs. Having one or more records for duplicated
  # variants will not change the number of patients, so the patient level
  # frequencies will not be chagned.
  #
  # Sample level frequencies: for each variant, check how many unique
  # Kids_First_Biospecimen_IDs are in primary or relapse independent
  # Kids_First_Biospecimen_IDs. Having one or more records for duplicated
  # variants will not change the number of Kids_First_Biospecimen_IDs, so the
  # sample level frequencies will not be chagned.
  sample_var_df <- ss_alt_df %>%
    select(Kids_First_Biospecimen_ID, Alt_ID) %>%
    distinct()
  
  patient_var_df <- sample_var_df %>%
    left_join(bpid_ss_htl_df, by = 'Kids_First_Biospecimen_ID') %>%
    group_by(Alt_ID) %>%
    summarise(Total_alterations = length(unique(Kids_First_Participant_ID))) %>%
    mutate(Patients_in_dataset = ss_n_patients) %>%
    mutate(Total_alterations_Over_Patients_in_dataset =
             paste(Total_alterations, Patients_in_dataset, sep = '/')) %>%
    mutate(Frequency_in_overall_dataset = num_to_pct_chr(
      Total_alterations / Patients_in_dataset))
  
  # td = tumor descriptor
  td_var_df <- sample_var_df %>%
    group_by(Alt_ID) %>%
    summarise(
      Total_primary_tumors_mutated =
        sum(unique(Kids_First_Biospecimen_ID) %in% ss_primary_kfbids),
      Total_relapse_tumors_mutated =
        sum(unique(Kids_First_Biospecimen_ID) %in% ss_relapse_kfbids)) %>%
    mutate(Primary_tumors_in_dataset = ss_n_primary_tumors,
           Relapse_tumors_in_dataset = ss_n_relapse_tumors) %>%
    mutate(Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset =
             paste(Total_primary_tumors_mutated,
                   Primary_tumors_in_dataset, sep = '/'),
           Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset =
             paste(Total_relapse_tumors_mutated, Relapse_tumors_in_dataset,
                   sep = '/')) %>%
    mutate(Frequency_in_primary_tumors = num_to_pct_chr(
      Total_primary_tumors_mutated / Primary_tumors_in_dataset)) %>%
    mutate(Frequency_in_relapse_tumors = num_to_pct_chr(
      Total_relapse_tumors_mutated / Relapse_tumors_in_dataset))
  
  output_var_df <- patient_var_df %>%
    left_join(td_var_df, by = 'Alt_ID') %>%
    arrange(desc(as.numeric(Total_alterations)))  %>%
    mutate(Disease = ss_cancer_group,
           Dataset = paste(ss_cohorts, collapse = '&')) %>%
    replace_na(list(
      Total_mutations_Over_Patients_in_dataset = '',
      Frequency_in_overall_dataset = '',
      Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = '',
      Frequency_in_primary_tumors = '',
      Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = '',
      Frequency_in_relapse_tumors = ''))
  
  return(output_var_df)
}
