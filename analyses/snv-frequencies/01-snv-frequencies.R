library(tidyverse)



# Function definitions ----------------------------------------------------
# Format numbers to percentage characters
#
# Adapted from @Richie Cotton's answer at
# https://stackoverflow.com/a/7146270/4638182
num_to_pct_chr <- function(x, digits = 2, format = "f", ...) {
  stopifnot(!is.null(x))
  if (length(x) == 0) {
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



# Collapse a list of refseq.protein character vectors from
# mygene.info query results
collapse_rp_lists <- function(xl) {
  fxl <- discard(xl, is.null)
  # assert all characters
  sapply(fxl, function(x) {
    stopifnot(is.character(x))
  })
  fxv <- unique(reduce(fxl, c, .init = character(0)))
  np_fxv <- fxv[str_detect(fxv, '^NP_')]
  c_np_fxv <- paste(np_fxv, collapse = ',')
  return(c_np_fxv)
}



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



# Generate overall/primary/relapse mutation frequency table for a set of samples
#
# Args:
# - maf_df: a MAF tibble. Must contain the following fields:
#   Kids_First_Biospecimen_ID and var_group_col. Kids_First_Biospecimen_IDs must
#   be a subset of the ones in overall_histology_df.
# - var_group_col: the column name in maf_df to group variants when computing
#   mutation frequencies. For example, if var_group_col is 'Variant_Unique_ID',
#   the mutation frequencies are at variant level; if var_group_col is 'Gene_ID'
#   the mutation frequences are at gene level.
# - overall_histology_df: the histology tibble that contains a set of samples
#   for computing mutation frequencies. Must contain the following fields:
#   Kids_First_Biospecimen_ID, Kids_First_Participant_ID. This is used for
#   computing the following columns: Total_mutations, Patients_in_dataset,
#   Total_mutations_Over_Patients_in_dataset and Frequency_in_overall_dataset.
# - primary_histology_df: the histology tibble that contains primary tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - relapse_histology_df: the histology tibble that contains relapse tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
#
# Returns a MAF tibble with additional frequency columns.
get_opr_mut_freq_tbl <- function(maf_df, var_group_col,
                                 overall_histology_df, primary_histology_df,
                                 relapse_histology_df) {
  # check input parameters
  stopifnot(is.character(var_group_col))
  stopifnot(identical(length(var_group_col), as.integer(1)))
  stopifnot(identical(sum(is.na(var_group_col)), as.integer(0)))
  stopifnot(var_group_col %in% colnames(maf_df))
  stopifnot('Kids_First_Biospecimen_ID' %in% colnames(maf_df))
  stopifnot('Kids_First_Biospecimen_ID' %in% colnames(overall_histology_df))
  stopifnot('Kids_First_Biospecimen_ID' %in% colnames(primary_histology_df))
  stopifnot('Kids_First_Biospecimen_ID' %in% colnames(relapse_histology_df))
  stopifnot('Kids_First_Participant_ID' %in% colnames(overall_histology_df))
  stopifnot('Kids_First_Participant_ID' %in% colnames(primary_histology_df))
  stopifnot('Kids_First_Participant_ID' %in% colnames(relapse_histology_df))
  # asser no NAs in Kids_First_Biospecimen_ID or var_group_col
  stopifnot(identical(
    sum(is.na(select_at(maf_df,
                        c('Kids_First_Biospecimen_ID', var_group_col)))),
    as.integer(0)
  ))
  stopifnot(identical(
    sum(is.na(select(overall_histology_df,
                     Kids_First_Biospecimen_ID,
                     Kids_First_Participant_ID))),
    as.integer(0)
  ))
  stopifnot(identical(
    sum(is.na(select(primary_histology_df,
                     Kids_First_Biospecimen_ID,
                     Kids_First_Participant_ID))),
    as.integer(0)
  ))
  stopifnot(identical(
    sum(is.na(select(relapse_histology_df,
                     Kids_First_Biospecimen_ID,
                     Kids_First_Participant_ID))),
    as.integer(0)
  ))
  # assert all samples in maf_df are in overall_histology_df
  stopifnot(all(maf_df$Kids_First_Biospecimen_ID %in%
                  overall_histology_df$Kids_First_Biospecimen_ID))

  # ss represents selected samples in overall_histology_df
  bpid_ss_htl_df <- overall_histology_df %>%
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

  # If one var_group_col value has more than one records in one
  # Kids_First_Biospecimen_ID, treat multiple records as a single mutation. The
  # following procedure is still correct when computing mutation frequencies.
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
  #
  # distinct() is necessary, because we only count one sample once, even if the
  # sample has multiple mutations that belong to the variant group
  sample_var_df <- maf_df %>%
    select_at(c('Kids_First_Biospecimen_ID', var_group_col)) %>%
    distinct()

  # Even if Kids_First_Biospecimen_ID has duplicates, the results will not be
  # affected, because only unique Kids_First_Biospecimen_ID and
  # Kids_First_Participant_ID are counted.
  patient_var_df <- sample_var_df %>%
    left_join(bpid_ss_htl_df, by = 'Kids_First_Biospecimen_ID') %>%
    group_by_at(var_group_col) %>%
    summarise(Total_mutations = length(unique(Kids_First_Participant_ID))) %>%
    mutate(Patients_in_dataset = ss_n_patients) %>%
    mutate(Total_mutations_Over_Patients_in_dataset =
             paste(Total_mutations, Patients_in_dataset, sep = '/')) %>%
    mutate(Frequency_in_overall_dataset = num_to_pct_chr(
      Total_mutations / Patients_in_dataset))

  # td = tumor descriptor
  td_var_df <- sample_var_df %>%
    group_by_at(var_group_col) %>%
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

  mut_freq_tbl <- inner_join(patient_var_df, td_var_df, by = var_group_col)
  # assertions to make sure the join is one-on-one
  stopifnot(identical(nrow(mut_freq_tbl), nrow(td_var_df)))
  stopifnot(identical(nrow(mut_freq_tbl), nrow(patient_var_df)))
  stopifnot(identical(
    nrow(mut_freq_tbl),
    length(unique(mut_freq_tbl[, var_group_col, drop = TRUE]))
  ))
  stopifnot(identical(sort(mut_freq_tbl[, var_group_col, drop = TRUE]),
                      sort(td_var_df[, var_group_col, drop = TRUE])))
  stopifnot(identical(sort(mut_freq_tbl[, var_group_col, drop = TRUE]),
                      sort(patient_var_df[, var_group_col, drop = TRUE])))
  stopifnot(identical(
    sort(colnames(mut_freq_tbl)),
    sort(c(var_group_col,
           'Total_mutations',
           'Patients_in_dataset',
           'Total_primary_tumors_mutated',
           'Total_relapse_tumors_mutated',
           'Primary_tumors_in_dataset',
           'Relapse_tumors_in_dataset',
           'Total_mutations_Over_Patients_in_dataset',
           'Frequency_in_overall_dataset',
           'Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset',
           'Frequency_in_primary_tumors',
           'Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset',
           'Frequency_in_relapse_tumors'))
  ))

  return(mut_freq_tbl)
}



# Return a single character value for representing a set of cohorts
#
# Args:
# - cohort_vec: a character vector of unique cohort values
#
# Returns a single character value for representing the cohort_vec.
get_cohort_set_value <- function(cohort_vec) {
  stopifnot(is.character(cohort_vec))
  stopifnot(!is.null(length(cohort_vec)))
  stopifnot(identical(length(cohort_vec), length(unique(cohort_vec))))
  stopifnot(identical(sum(is.na(cohort_vec)), as.integer(0)))
  stopifnot(length(cohort_vec) >= 1)

  if (identical(length(cohort_vec), as.integer(1))) {
    return(cohort_vec)
  } else {
    return('all_cohorts')
  }
}



# Convert cohort, cancer_group to PedcBio PedOT case_set_id
#
# Args:
# - ss_cancer_group: a single character value of cancer_group
# - cohort_set: a character vector of a set of unique cohort values
#
# Returns a single character value of case_set_id
#
# Note: the conversion procedure is by convention discussed at
# https://github.com/PediatricOpenTargets/ticket-tracker/issues/93
# #issuecomment-877397386
#
# Convertion rules:
# 1. cancer_group
#   i. retain alphanum, -, _
#   ii. other chars to _
#   iii. to lower acse
# 2. cohort_set
#   - if length 1:
#     i. retain alphanum, -, _
#     ii. other chars to _
#     iii. to lower acse
#     iiii. concatente to cancer_group by "_"
#   - if length > 1: drop
# 3. prepend ped_opentargets_2021_
get_pcb_pot_csi <- function(ss_cancer_group, cohort_set) {
  stopifnot(!is.null(ss_cancer_group))
  stopifnot(is.character(ss_cancer_group))
  stopifnot(identical(length(ss_cancer_group), as.integer(1)))
  stopifnot(identical(sum(is.na(ss_cancer_group)), as.integer(0)))

  # cohort_set is necessary to check whether one or more cohorts are used
  # cohort_set cannot have 0 length
  stopifnot(is.character(cohort_set))
  stopifnot(!is.null(length(cohort_set)))
  stopifnot(identical(length(cohort_set), length(unique(cohort_set))))
  stopifnot(identical(sum(is.na(cohort_set)), as.integer(0)))
  stopifnot(length(cohort_set) >= 1)

  cancer_group_stable <- tolower(
    gsub("[^-_a-zA-Z0-9]", "_", ss_cancer_group))

  if (length(cohort_set) == 1) {
    # only one cohort
    # add cohort value to case_set_id
    cohort_stable <- tolower(gsub("[^-_a-zA-Z0-9]", "_", cohort_set))

    case_set_id <- paste0('ped_opentargets_2021_', cohort_stable, '_',
                          cancer_group_stable)
  } else {
    # more than one cohort
    # omit cohort value in case_set_id
    case_set_id <- paste0('ped_opentargets_2021_', cancer_group_stable)
  }

  return(case_set_id)
}



# Convert gene symbol, cohort, cancer_group, and plot type to PedcBio PedOT plot
# URL
#
# Args:
# - gene_symbol_vec: a character vector of gene symbols
# - ss_case_set_id: a single character value of PedcBio PedOT case_set_id
# - plot_type: a single character value following
#   https://pedcbioportal.kidsfirstdrc.org/results/ in the URL, e.g. 'oncoprint'
#   and 'mutations'.
#
# Returns a character vector of oncoprint plot URL
get_pcb_pot_plot_url <- function(gene_symbol_vec, ss_case_set_id, plot_type) {
  stopifnot(!is.null(gene_symbol_vec))
  stopifnot(is.character(gene_symbol_vec))
  stopifnot(identical(sum(is.na(gene_symbol_vec)), as.integer(0)))

  stopifnot(!is.null(ss_case_set_id))
  stopifnot(is.character(ss_case_set_id))
  stopifnot(identical(sum(is.na(ss_case_set_id)), as.integer(0)))
  stopifnot(identical(length(ss_case_set_id), as.integer(1)))

  stopifnot(!is.null(plot_type))
  stopifnot(is.character(plot_type))
  stopifnot(identical(sum(is.na(plot_type)), as.integer(0)))
  stopifnot(identical(length(plot_type), as.integer(1)))

  if (length(gene_symbol_vec) == 0) {
    return(character(0))
  }

  plot_url_vec <- paste0(
    'https://pedcbioportal.kidsfirstdrc.org/results/', plot_type,
    '?cancer_study_list=ped_opentargets_2021&case_set_id=',
    ss_case_set_id,
    '&Action=Submit&gene_list=', gene_symbol_vec)

  return(plot_url_vec)
}



# Add PedcBio PedOT oncoprint and mutations plot URLs to the mutation frequency
# table
#
# Args:
# - mut_freq_tbl: a mutation frequency tibble. Must contain Gene_symbol column.
# - ss_cancer_group: a character value of the cancer group to compute for.
# - ss_cohorts: a vector of character values of the cohorts to compute for.
# - valid_url_case_set_ids: a vector of character values that are valid
#   case_set_ids for generating PedcBio PedOT oncoprint and mutations plot URLs.
#
# Returns a mutation frequency tibble with additional
# PedcBio_PedOT_oncoprint_plot_URL and PedcBio_PedOT_mutations_plot_URL columns
add_cg_ch_pedcbio_pedot_plot_urls <- function(mut_freq_tbl,
                                              ss_cancer_group,
                                              ss_cohorts,
                                              valid_url_case_set_ids) {
  # check input parameters
  stopifnot(is.character(ss_cancer_group))
  stopifnot(identical(length(ss_cancer_group), as.integer(1)))
  stopifnot(identical(sum(is.na(ss_cancer_group)), as.integer(0)))

  stopifnot(is.character(ss_cohorts))
  stopifnot(!is.null(length(ss_cohorts)))
  stopifnot(identical(length(ss_cohorts), length(unique(ss_cohorts))))
  stopifnot(identical(sum(is.na(ss_cohorts)), as.integer(0)))
  stopifnot(length(ss_cohorts) >= 1)

  stopifnot(is.character(valid_url_case_set_ids))

  # a = annotated
  ss_case_set_id <- get_pcb_pot_csi(ss_cancer_group, ss_cohorts)
  if (ss_case_set_id %in% valid_url_case_set_ids) {
    a_mut_freq_tbl <- mut_freq_tbl %>%
      mutate(PedcBio_PedOT_oncoprint_plot_URL = get_pcb_pot_plot_url(
        Gene_symbol, ss_case_set_id, 'oncoprint')) %>%
      mutate(PedcBio_PedOT_mutations_plot_URL = get_pcb_pot_plot_url(
        Gene_symbol, ss_case_set_id, 'mutations'))
  } else {
    a_mut_freq_tbl <- mut_freq_tbl %>%
      mutate(PedcBio_PedOT_oncoprint_plot_URL = '') %>%
      mutate(PedcBio_PedOT_mutations_plot_URL = '')
  }

  return(a_mut_freq_tbl)
}



# Generate variant-level mutation frequency table for a (cancer_group, cohort)
#
# Args:
# - maf_df: a MAF tibble. Must contain the following fields:
#   Kids_First_Biospecimen_ID,
#   Variant_ID,
#   Hugo_Symbol,
#   Gene_full_name,
#   dbSNP_RS,
#   IMPACT,
#   SIFT,
#   PolyPhen,
#   Variant_Classification,
#   Variant_Type,
#   Protein_RefSeq_ID,
#   Gene,
#   ENSP,
#   HGVSp_Short,
#   HotSpotAllele.
# - overall_histology_df: the histology tibble that contains all samples. Must
#   contain the following fields: Kids_First_Biospecimen_ID,
#   Kids_First_Participant_ID, cancer_group, and cohort. This
#   is used for computing the following columns: Total_mutations,
#   Patients_in_dataset, Total_mutations_Over_Patients_in_dataset and
#   Frequency_in_overall_dataset.
# - primary_histology_df: the histology tibble that contains primary tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - relapse_histology_df: the histology tibble that contains relapse tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - ss_cancer_group: a character value of the cancer group to compute for.
# - ss_cohorts: a vector of character values of the cohorts to compute for.
# - valid_url_case_set_ids: a vector of character values that are valid
#   case_set_ids for generating PedcBio PedOT oncoprint and mutations plot URLs.
#
# Returns a MAF tibble with additional frequency columns.
get_cg_ch_var_level_mut_freq_tbl <- function(maf_df, overall_histology_df,
                                             primary_histology_df,
                                             relapse_histology_df,
                                             ss_cancer_group, ss_cohorts,
                                             valid_url_case_set_ids) {
  # check input parameters
  stopifnot(is.character(ss_cancer_group))
  stopifnot(identical(length(ss_cancer_group), as.integer(1)))
  stopifnot(identical(sum(is.na(ss_cancer_group)), as.integer(0)))

  stopifnot(is.character(ss_cohorts))
  stopifnot(!is.null(length(ss_cohorts)))
  stopifnot(identical(length(ss_cohorts), length(unique(ss_cohorts))))
  stopifnot(identical(sum(is.na(ss_cohorts)), as.integer(0)))
  stopifnot(length(ss_cohorts) >= 1)

  # ss = subset
  ss_htl_df <- overall_histology_df %>%
    filter(cancer_group == ss_cancer_group,
           cohort %in% ss_cohorts)

  ss_maf_df <- maf_df %>%
    filter(Kids_First_Biospecimen_ID %in% ss_htl_df$Kids_First_Biospecimen_ID)
  # Need to subset overall_histology_df and maf_df, because they are not subset
  # in get_opr_mut_freq_tbl().
  #
  # No need to subset primary_histology_df and relapse_histology_df, because
  # they are subset in get_opr_mut_freq_tbl().
  #
  # ss_mut_freq_df Variant_ID is asserted all unique in get_opr_mut_freq_tbl
  ss_mut_freq_df <- get_opr_mut_freq_tbl(ss_maf_df, 'Variant_ID',
                                         ss_htl_df, primary_histology_df,
                                         relapse_histology_df)

  # If one Variant_ID has more than 1 values in the summarised fields, add
  # code to handle duplicates.
  # In case we need mRNA_RefSeq_ID in the future, add to the summarise call.
  #            mRNA_RefSeq_ID = unique(RefSeq),
  # Allow one Variant_ID mapped to multiple ENSP IDs, in order to be consistent
  # with the gene-level table.
  output_var_df <- ss_maf_df %>%
    group_by(Variant_ID) %>%
    summarise(Gene_symbol = unique(Hugo_Symbol),
              Gene_full_name = unique(Gene_full_name),
              dbSNP_ID = unique(dbSNP_RS),
              VEP_impact = unique(IMPACT),
              SIFT_impact = unique(SIFT),
              PolyPhen_impact = unique(PolyPhen),
              Variant_classification = unique(Variant_Classification),
              Variant_type = unique(Variant_Type),
              Protein_RefSeq_ID = unique(Protein_RefSeq_ID),
              Gene_Ensembl_ID = unique(Gene),
              Protein_Ensembl_ID = paste(
                discard(unique(ENSP), is.na), collapse = ','),
              Protein_change = unique(HGVSp_Short),
              HotSpotAllele = unique(HotSpotAllele)) %>%
    left_join(ss_mut_freq_df, by = 'Variant_ID') %>%
    mutate(Disease = ss_cancer_group,
           Dataset = get_cohort_set_value(ss_cohorts),
           HotSpot = if_else(HotSpotAllele == 1, true = 'Y', false = 'N')) %>%
    arrange(desc(as.numeric(Total_mutations))) %>%
    select(Gene_symbol, Dataset, Disease, Variant_ID, dbSNP_ID,
           VEP_impact, SIFT_impact, PolyPhen_impact, Variant_classification,
           Variant_type, Gene_full_name, Protein_RefSeq_ID,
           Gene_Ensembl_ID, Protein_Ensembl_ID, Protein_change,
           Total_mutations_Over_Patients_in_dataset,
           Frequency_in_overall_dataset,
           Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
           Frequency_in_primary_tumors,
           Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
           Frequency_in_relapse_tumors, HotSpot) %>%
    mutate_all(function(x) replace_na(x, replace = ''))

  output_var_df <- add_cg_ch_pedcbio_pedot_plot_urls(
    output_var_df, ss_cancer_group, ss_cohorts, valid_url_case_set_ids)
  stopifnot(identical(sum(is.na(output_var_df)), as.integer(0)))

  return(output_var_df)
}



# Generate gene-level mutation frequency table for a (cancer_group, cohort)
#
# Args:
# - maf_df: a MAF tibble. Must contain the following fields:
#   Kids_First_Biospecimen_ID,
#   Hugo_Symbol,
#   Gene_full_name,
#   Protein_RefSeq_ID,
#   Gene,
#   ENSP.
# - overall_histology_df: the histology tibble that contains all samples. Must
#   contain the following fields: Kids_First_Biospecimen_ID,
#   Kids_First_Participant_ID, cancer_group, and cohort. This
#   is used for computing the following columns: Total_mutations,
#   Patients_in_dataset, Total_mutations_Over_Patients_in_dataset and
#   Frequency_in_overall_dataset.
# - primary_histology_df: the histology tibble that contains primary tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - relapse_histology_df: the histology tibble that contains relapse tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - ss_cancer_group: a single character value of the cancer group to compute
#   for.
# - ss_cohorts: a vector of character values of the cohorts to compute for.
# - valid_url_case_set_ids: a vector of character values that are valid
#   case_set_ids for generating PedcBio PedOT oncoprint and mutations plot URLs.
#
# Returns a MAF tibble with additional frequency columns.
get_cg_ch_gene_level_mut_freq_tbl <- function(maf_df, overall_histology_df,
                                              primary_histology_df,
                                              relapse_histology_df,
                                              ss_cancer_group, ss_cohorts,
                                              valid_url_case_set_ids) {
  # check input parameters
  stopifnot(is.character(ss_cancer_group))
  stopifnot(identical(length(ss_cancer_group), as.integer(1)))
  stopifnot(identical(sum(is.na(ss_cancer_group)), as.integer(0)))

  stopifnot(is.character(ss_cohorts))
  stopifnot(!is.null(length(ss_cohorts)))
  stopifnot(identical(length(ss_cohorts), length(unique(ss_cohorts))))
  stopifnot(identical(sum(is.na(ss_cohorts)), as.integer(0)))
  stopifnot(length(ss_cohorts) >= 1)

  # ss = subset
  ss_htl_df <- overall_histology_df %>%
    filter(cancer_group == ss_cancer_group,
           cohort %in% ss_cohorts)

  ss_maf_df <- maf_df %>%
    filter(Kids_First_Biospecimen_ID %in% ss_htl_df$Kids_First_Biospecimen_ID)
  # Need to subset overall_histology_df and maf_df, because they are not subset
  # in get_opr_mut_freq_tbl().
  #
  # No need to subset primary_histology_df and relapse_histology_df, because
  # they are subset in get_opr_mut_freq_tbl().
  #
  # Gene column is Gene Ensembl ENSG ID
  #
  # ss_mut_freq_df Variant_ID is asserted all unique in get_opr_mut_freq_tbl
  ss_mut_freq_df <- get_opr_mut_freq_tbl(ss_maf_df, 'Gene',
                                         ss_htl_df, primary_histology_df,
                                         relapse_histology_df)

  # If one Gene has more than 1 values in the summarised fields, add
  # code to handle duplicates.
  # In case we need mRNA_RefSeq_ID in the future, add to the summarise call.
  #            mRNA_RefSeq_ID = unique(RefSeq),
  output_var_df <- ss_maf_df %>%
    group_by(Gene) %>%
    summarise(Gene_symbol = unique(Hugo_Symbol),
              Gene_full_name = unique(Gene_full_name),
              Protein_RefSeq_ID = unique(Protein_RefSeq_ID),
              Protein_Ensembl_ID = paste(
                discard(unique(ENSP), is.na), collapse = ',')) %>%
    left_join(ss_mut_freq_df, by = 'Gene') %>%
    rename(Gene_Ensembl_ID = Gene) %>%
    mutate(Disease = ss_cancer_group,
           Dataset = get_cohort_set_value(ss_cohorts)) %>%
    arrange(desc(as.numeric(Total_mutations))) %>%
    select(Gene_symbol, Dataset, Disease, Gene_full_name, Protein_RefSeq_ID,
           Gene_Ensembl_ID, Protein_Ensembl_ID,
           Total_mutations_Over_Patients_in_dataset,
           Frequency_in_overall_dataset,
           Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
           Frequency_in_primary_tumors,
           Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
           Frequency_in_relapse_tumors) %>%
    mutate_all(function(x) replace_na(x, replace = ''))

  output_var_df <- add_cg_ch_pedcbio_pedot_plot_urls(
    output_var_df, ss_cancer_group, ss_cohorts, valid_url_case_set_ids)
  stopifnot(identical(sum(is.na(output_var_df)), as.integer(0)))

  return(output_var_df)
}



# Create output dir ------------------------------------------------------------
tables_dir <- 'results'

if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)
}



# Read data --------------------------------------------------------------------
message('Read data...')
htl_df <- read_tsv('../../data/histologies.tsv', guess_max = 100000,
                   col_types = cols(.default = col_guess()))
# assert no Kids_First_Biospecimen_ID or Kids_First_Participant_ID is NA
stopifnot(identical(
  sum(is.na(select(htl_df, Kids_First_Biospecimen_ID,
                   Kids_First_Participant_ID))),
  as.integer(0)))

maf_df <- read_tsv(
  '../../data/snv-consensus-plus-hotspots.maf.tsv.gz', comment = '#',
  col_types = cols(
    .default = col_guess(),
    CLIN_SIG = col_character(),
    PUBMED = col_character()))
# assert all NCBI_Build values are GRCh38
stopifnot(all(maf_df$NCBI_Build == 'GRCh38'))
# assert all records have tumor sample barcode
stopifnot(identical(sum(is.na(maf_df$Tumor_Sample_Barcode)), as.integer(0)))

# primary independent sample data frame
primary_indp_sdf <- read_tsv(
  file.path('..', 'independent-samples', 'results',
            'independent-specimens.wgswxspanel.primary.eachcohort.tsv'),
  col_types = cols(
    .default = col_guess()))

# relapse independent samples
relapse_indp_sdf <- read_tsv(
  file.path('..', 'independent-samples', 'results',
            'independent-specimens.wgswxspanel.relapse.eachcohort.tsv'),
  col_types = cols(
    .default = col_guess()))

# efo cancer_group mappings
efo_mondo_cg_df <- read_tsv('../../data/efo-mondo-map.tsv',
                            col_types = cols(.default = col_guess())) %>%
  distinct()
# efo_mondo_cg_df$cancer_group[
#   !efo_mondo_cg_df$cancer_group %in% htl_df$cancer_group]

# assert all cancer_groups are not NA
stopifnot(identical(sum(is.na(efo_mondo_cg_df$cancer_group)), as.integer(0)))
# assert all cancer_groups are unique.
# result SNV table is left joined by cancer_groups
stopifnot(identical(length(unique(efo_mondo_cg_df$cancer_group)),
                    nrow(efo_mondo_cg_df)))

# ensg hugo rmtl mappings
ensg_hugo_rmtl_df <- read_tsv('../../data/ensg-hugo-rmtl-v1-mapping.tsv',
                              col_types = cols(.default = col_guess())) %>%
  distinct()
# assert all ensg_ids and gene_symbols are not NA
stopifnot(identical(sum(is.na(ensg_hugo_rmtl_df$ensg_id)), as.integer(0)))
stopifnot(identical(sum(is.na(ensg_hugo_rmtl_df$gene_symbol)), as.integer(0)))
# assert all ensg_id are unique
# result SNV table is left joined by ensg_id
stopifnot(identical(length(unique(ensg_hugo_rmtl_df$ensg_id)),
                    nrow(ensg_hugo_rmtl_df)))

# gene symbol and gene type mapping
gsb_gtype_df <- read_tsv('../fusion_filtering/references/genelistreference.txt',
                         col_types = cols(.default = col_guess()))
# assert no NA in gsb_gtype_df
stopifnot(identical(sum(is.na(gsb_gtype_df)), as.integer(0)))

# pcb = PedcBioPortal
# pot = Pediatric Open Targets
pcb_pot_case_set_list <- jsonlite::read_json(
  'input/ped_opentargets_2021_pedcbio_case_set_ids.json')

pcb_pot_case_set_id_vec <- vapply(
  pcb_pot_case_set_list,
  function(x) {
    stopifnot(!is.null(x))
    stopifnot(!is.null(x$sampleListId))
    stopifnot(is.character(x$sampleListId))
    stopifnot(identical(length(x$sampleListId), as.integer(1)))
    stopifnot(identical(sum(is.na(x$sampleListId)), as.integer(0)))
    return(x$sampleListId)
  },
  FUN.VALUE = character(1)
)
stopifnot(identical(length(pcb_pot_case_set_id_vec),
                    length(unique(pcb_pot_case_set_id_vec))))

# gene symbol to OncoKB cancer gene, oncogene and tumor suppressor gene mappings
oncokb_cancer_gene_df <- read_tsv('input/oncokb_cancer_gene_list.tsv',
                                  col_types = cols(.default = col_guess()))
# assert no NA in gene symbols
stopifnot(identical(sum(is.na(select(oncokb_cancer_gene_df, `Hugo Symbol`))),
                    as.integer(0)))
# assert all symbols are unique
stopifnot(identical(nrow(oncokb_cancer_gene_df),
                    length(unique(oncokb_cancer_gene_df$`Hugo Symbol`))))


# Subset tumor samples and used columns in MAF table ---------------------------
tumor_kfbids <- htl_df %>%
  filter(sample_type == 'Tumor', !is.na(cancer_group),
         !is.na(cohort)) %>%
  pull(Kids_First_Biospecimen_ID)

# It is ok to subset non-synonymous variants before computing frequencies,
# because the total numbers of patients, primary samples, and relapse samples
# are determined by histology dataframes.
#
# If a sample has only synonymous variants, it will still be counted in the
# total numbers, but it will not be counted in any mutated number of
# non-synonymous mutation.
maf_df <- maf_df %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del",
                                       "Frame_Shift_Ins",
                                       "Splice_Site",
                                       "Nonsense_Mutation",
                                       "Nonstop_Mutation",
                                       "In_Frame_Del",
                                       "In_Frame_Ins",
                                       "Missense_Mutation",
                                       "Fusion",
                                       "Multi_Hit",
                                       "Multi_Hit_Fusion",
                                       "Hom_Deletion",
                                       "Hem_Deletion",
                                       "Amp",
                                       "Del",
                                       "Translation_Start_Site")) %>%
  filter(Tumor_Sample_Barcode %in% tumor_kfbids) %>%
  mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  mutate(Variant_ID = paste(Chromosome, Start_Position, Reference_Allele,
                            Tumor_Seq_Allele2, sep = '_')) %>%
  select(Kids_First_Biospecimen_ID, Variant_ID,
         Hugo_Symbol, dbSNP_RS, IMPACT, SIFT, PolyPhen, Variant_Classification,
         Variant_Type, RefSeq, Gene, ENSP, HGVSp_Short, HotSpotAllele)

rm(tumor_kfbids)



# Subset independent samples in histology table --------------------------------
message('Prepare data...')
# subset histology dataframe to have only samples that are in SNV MAF
maf_sample_htl_df <- htl_df %>%
  filter(Kids_First_Biospecimen_ID %in% maf_df$Kids_First_Biospecimen_ID)
# assert all samples are unique
stopifnot(identical(
  nrow(maf_sample_htl_df),
  length(unique(maf_sample_htl_df$Kids_First_Biospecimen_ID))
))
# assert samples are the same
stopifnot(identical(sort(unique(maf_df$Kids_First_Biospecimen_ID)),
                    sort(maf_sample_htl_df$Kids_First_Biospecimen_ID)))

# td = tumor descriptor
td_htl_dfs <- list(
  primary_htl_df = maf_sample_htl_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             primary_indp_sdf$Kids_First_Biospecimen_ID),

  relapse_htl_df = maf_sample_htl_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             relapse_indp_sdf$Kids_First_Biospecimen_ID),

  overall_htl_df = maf_sample_htl_df
)



# Add additional annotations ---------------------------------------------------
message('Retrieve Gene_full_name and Protein_RefSeq_ID from mygene.info...')
maf_ens_gids <- unique(maf_df$Gene)
maf_ens_gids <- maf_ens_gids[!is.na(maf_ens_gids)]

mg_qres_list <- mygene::queryMany(
  maf_ens_gids, scopes = 'ensembl.gene', fields = c('refseq', 'name'),
  species = 'human', returnall = TRUE, return.as = 'DataFrame')

mg_qres_df <- as_tibble(
  mg_qres_list$response[, c('query', 'notfound', 'name', 'refseq.protein')]) %>%
  replace_na(list(notfound = FALSE)) %>%
  filter(!notfound) %>%
  group_by(query) %>%
  summarise(name = paste(unique(name), collapse = ', '),
            refseq_protein = collapse_rp_lists(refseq.protein)) %>%
  rename(Gene = query, Gene_full_name = name,
         Protein_RefSeq_ID = refseq_protein)

# add additional fields to MAF df
# mg_qres_df$Gene values are unique, because they are the group_by column
maf_df <- maf_df %>%
  left_join(mg_qres_df, by = 'Gene')



# Compute mutation frequencies -------------------------------------------------
message('Compute mutation frequencies...')
cancer_group_cohort_summary_df <- get_cg_cs_tbl(td_htl_dfs$overall_htl_df)

# nf = n_samples filtered
nf_cancer_group_cohort_summary_df <- cancer_group_cohort_summary_df %>%
  filter(n_samples >= 5)

mut_freq_tbl_list <- lapply(
  seq_len(nrow(nf_cancer_group_cohort_summary_df)),
  function(i) {
    cgcs_row <- nf_cancer_group_cohort_summary_df[i, ]
    stopifnot(identical(nrow(cgcs_row), as.integer(1)))
    stopifnot(identical(length(cgcs_row$cohort_list), as.integer(1)))
    stopifnot(is.character(cgcs_row$cohort_list[[1]]))

    c_cancer_group <- cgcs_row$cancer_group
    c_cohorts <- cgcs_row$cohort_list[[1]]
    stopifnot(identical(paste(c_cohorts, collapse = '&'), cgcs_row$cohort))
    message(paste(c_cancer_group, cgcs_row$cohort))

    var_level_tbl <- get_cg_ch_var_level_mut_freq_tbl(
      maf_df, td_htl_dfs$overall_htl_df, td_htl_dfs$primary_htl_df,
      td_htl_dfs$relapse_htl_df, c_cancer_group, c_cohorts,
      pcb_pot_case_set_id_vec)

    gene_level_tbl <- get_cg_ch_gene_level_mut_freq_tbl(
      maf_df, td_htl_dfs$overall_htl_df, td_htl_dfs$primary_htl_df,
      td_htl_dfs$relapse_htl_df, c_cancer_group, c_cohorts,
      pcb_pot_case_set_id_vec)

    res_list <- list(var_level_tbl = var_level_tbl,
                     gene_level_tbl = gene_level_tbl)
    return(res_list)
  }
)

stopifnot(identical(length(mut_freq_tbl_list),
                    nrow(nf_cancer_group_cohort_summary_df)))

var_level_mut_freq_tbl <- distinct(bind_rows(
  lapply(mut_freq_tbl_list, function(x) x$var_level_tbl)
))
# assert all rows are unique
stopifnot(identical(
  nrow(var_level_mut_freq_tbl),
  sum(sapply(mut_freq_tbl_list, function(x) nrow(x$var_level_tbl)))
))

gene_level_mut_freq_tbl <- distinct(bind_rows(
  lapply(mut_freq_tbl_list, function(x) x$gene_level_tbl)
))
# assert all rows are unique
stopifnot(identical(
  nrow(gene_level_mut_freq_tbl),
  sum(sapply(mut_freq_tbl_list, function(x) nrow(x$gene_level_tbl)))
))



# Add annotations to the output table ------------------------------------------
# Add EFO and MONDO
ann_efo_mondo_cg_df <- efo_mondo_cg_df %>%
  rename(Disease = cancer_group, EFO = efo_code, MONDO = mondo_code)

# unieuqness of efo_mondo_cg_df$cancer_group is asserted after reading
var_level_mut_freq_tbl <- var_level_mut_freq_tbl %>%
  left_join(ann_efo_mondo_cg_df, by = 'Disease') %>%
  replace_na(list(EFO = '', MONDO = ''))
stopifnot(identical(sum(is.na(var_level_mut_freq_tbl)), as.integer(0)))

gene_level_mut_freq_tbl <- gene_level_mut_freq_tbl %>%
  left_join(ann_efo_mondo_cg_df, by = 'Disease') %>%
  replace_na(list(EFO = '', MONDO = ''))
stopifnot(identical(sum(is.na(gene_level_mut_freq_tbl)), as.integer(0)))


# Add RMTL
# asert all rmtl NAs have version NAs, vice versa
stopifnot(identical(is.na(ensg_hugo_rmtl_df$rmtl),
                    is.na(ensg_hugo_rmtl_df$version)))
ann_ensg_hugo_rmtl_df <- ensg_hugo_rmtl_df %>%
  select(ensg_id, rmtl, version) %>%
  filter(!is.na(rmtl), !is.na(version)) %>%
  mutate(RMTL = paste0(rmtl, ' (', version, ')')) %>%
  select(ensg_id, RMTL) %>%
  rename(Gene_Ensembl_ID = ensg_id)

# unieuqness of ensg_hugo_rmtl_df$ensg_id is asserted after reading
var_level_mut_freq_tbl <- var_level_mut_freq_tbl %>%
  left_join(ann_ensg_hugo_rmtl_df, by = 'Gene_Ensembl_ID') %>%
  replace_na(list(RMTL = ''))
stopifnot(identical(sum(is.na(var_level_mut_freq_tbl)), as.integer(0)))

gene_level_mut_freq_tbl <- gene_level_mut_freq_tbl %>%
  left_join(ann_ensg_hugo_rmtl_df, by = 'Gene_Ensembl_ID') %>%
  replace_na(list(RMTL = ''))
stopifnot(identical(sum(is.na(gene_level_mut_freq_tbl)), as.integer(0)))


# Add gene type
ann_gsb_gtype_df <- gsb_gtype_df %>%
  group_by(Gene_Symbol) %>%
  summarise(Gene_type = paste(sort(unique(discard(type, is.na))),
                              collapse = ',')) %>%
  rename(Gene_symbol = Gene_Symbol)

# ann_gsb_gtype_df$Gene_type values are unique, because they are the group_by
# column
#
# Only add gene type column to gene-level table

gene_level_mut_freq_tbl <- gene_level_mut_freq_tbl %>%
  left_join(ann_gsb_gtype_df, by = 'Gene_symbol') %>%
  replace_na(list(Gene_type = ''))
stopifnot(identical(sum(is.na(gene_level_mut_freq_tbl)), as.integer(0)))

# Add OncoKB cancer gene and oncogene/TSG (tumor suppressor gene)
# all genes in input/oncokb_cancer_gene_list.tsv are OncoKB cancer genes
oncokb_cancer_gene_ann_df <- oncokb_cancer_gene_df %>%
  select(`Hugo Symbol`, `Is Oncogene`, `Is Tumor Suppressor Gene`) %>%
  rename(Gene_symbol = `Hugo Symbol`,
         is_onco = `Is Oncogene`,
         is_tsg = `Is Tumor Suppressor Gene`) %>%
  mutate(OncoKB_cancer_gene = 'Y',
         OncoKB_oncogene_TSG = case_when(
           is_onco == 'Yes' & is_tsg == 'Yes' ~ 'Oncogene,TumorSuppressorGene',
           is_onco == 'Yes' ~ 'Oncogene',
           is_tsg == 'Yes' ~ 'TumorSuppressorGene',
           TRUE ~ '')) %>%
  select(Gene_symbol, OncoKB_cancer_gene, OncoKB_oncogene_TSG)
# unieuqness of oncokb_cancer_gene_df$`Hugo Symbol` is asserted after reading
var_level_mut_freq_tbl <- var_level_mut_freq_tbl %>%
  left_join(oncokb_cancer_gene_ann_df, by = 'Gene_symbol') %>%
  replace_na(list(OncoKB_cancer_gene = 'N', OncoKB_oncogene_TSG = ''))
stopifnot(identical(sum(is.na(var_level_mut_freq_tbl)), as.integer(0)))

gene_level_mut_freq_tbl <- gene_level_mut_freq_tbl %>%
  left_join(oncokb_cancer_gene_ann_df, by = 'Gene_symbol') %>%
  replace_na(list(OncoKB_cancer_gene = 'N', OncoKB_oncogene_TSG = ''))
stopifnot(identical(sum(is.na(gene_level_mut_freq_tbl)), as.integer(0)))




# Output tsv and JSON -----------------------------------------------------
# reorder for output
var_level_mut_freq_tbl <- var_level_mut_freq_tbl %>%
  select(Gene_symbol, RMTL, Dataset, Disease, EFO, MONDO, Variant_ID, dbSNP_ID,
         VEP_impact, SIFT_impact, PolyPhen_impact, Variant_classification,
         Variant_type, Gene_full_name, Protein_RefSeq_ID,
         Gene_Ensembl_ID, Protein_Ensembl_ID, Protein_change,
         Total_mutations_Over_Patients_in_dataset,
         Frequency_in_overall_dataset,
         Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
         Frequency_in_primary_tumors,
         Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
         Frequency_in_relapse_tumors,
         HotSpot, OncoKB_cancer_gene, OncoKB_oncogene_TSG,
         PedcBio_PedOT_oncoprint_plot_URL, PedcBio_PedOT_mutations_plot_URL) %>%
  rename(Variant_ID_hg38 = Variant_ID)

gene_level_mut_freq_tbl <- gene_level_mut_freq_tbl %>%
  select(Gene_symbol, RMTL, Dataset, Disease, EFO, MONDO,
         Gene_full_name, Gene_type, Protein_RefSeq_ID,
         Gene_Ensembl_ID, Protein_Ensembl_ID,
         Total_mutations_Over_Patients_in_dataset,
         Frequency_in_overall_dataset,
         Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
         Frequency_in_primary_tumors,
         Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
         Frequency_in_relapse_tumors,
         OncoKB_cancer_gene, OncoKB_oncogene_TSG,
         PedcBio_PedOT_oncoprint_plot_URL, PedcBio_PedOT_mutations_plot_URL)

write_tsv(
  var_level_mut_freq_tbl,
  file.path(tables_dir, 'variant-level-snv-consensus-annotated-mut-freq.tsv'))

jsonlite::write_json(
  var_level_mut_freq_tbl,
  file.path(tables_dir, 'variant-level-snv-consensus-annotated-mut-freq.json'))

write_tsv(
  gene_level_mut_freq_tbl,
  file.path(tables_dir, 'gene-level-snv-consensus-annotated-mut-freq.tsv'))

jsonlite::write_json(
  gene_level_mut_freq_tbl,
  file.path(tables_dir, 'gene-level-snv-consensus-annotated-mut-freq.json'))


message('Done running 01-snv-frequencies.R.')
