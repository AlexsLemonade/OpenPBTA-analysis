suppressPackageStartupMessages(library(tidyverse))


# Generate means, standard deviations, z-scores, and ranks within each group.
#
# Args:
# - exp_df: (n_genes, n_samples) expression level data frame
# - groups: character vector of character vector of length n_samples
#
# Returns a list of (n_genes, n_groups) summary statistics tables.
get_expression_summary_stats <- function(exp_df, groups) {
  # unique groups to check that the computing steps do not modify the groups.
  check_groups <- sort(unique(groups))
  # gene symbols to check that the computing steps do not modify
  # the rownames of the exp_df.
  check_gids <- rownames(exp_df)

  # set check.names = FALSE and check.rows = FALSE to avoid R
  # changes the rownames or colnames implicitly
  c_exp_df <- data.frame(t(exp_df), check.names = FALSE,
                         check.rows = FALSE)
  c_exp_df$sample_group <- groups

  res_list <- list()


  # computed group means
  print('Compute means...')
  cg_mean_exp_df <- c_exp_df %>%
    group_by(sample_group) %>%
    summarise_all(mean) %>%
    column_to_rownames('sample_group')

  cg_mean_exp_out_df <- data.frame(
    t(cg_mean_exp_df), check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same
  stopifnot(identical(rownames(cg_mean_exp_out_df), check_gids))
  # assert groups are the same. assumes the first column is gene
  stopifnot(identical(sort(colnames(cg_mean_exp_out_df)),
                      check_groups))
  res_list$mean_df <- cg_mean_exp_out_df


  # compute group standard deviations
  print('Compute standard deviations...')
  cg_sd_exp_df <- c_exp_df %>%
    group_by(sample_group) %>%
    summarise_all(sd) %>%
    column_to_rownames('sample_group')

  cg_sd_exp_out_df <- data.frame(
    t(cg_sd_exp_df), check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same
  stopifnot(identical(rownames(cg_sd_exp_out_df), check_gids))
  # assert groups are the same. assumes the first column is gene
  stopifnot(identical(sort(colnames(cg_sd_exp_out_df)),
                      check_groups))
  res_list$sd_df <- cg_sd_exp_out_df


  # compute z-scores
  print('Compute z-scores...')
  # cg_mean_exp_df is (n_groups, n_genes)
  cg_mean_exp_mat <- as.matrix(cg_mean_exp_df)
  # assert gene symbols are the same
  stopifnot(identical(colnames(cg_mean_exp_mat), check_gids))
  # assert cancer groups are the same
  stopifnot(identical(sort(rownames(cg_mean_exp_mat)),
                      check_groups))
  # group wise means
  cg_mean_cgw_means <- rowMeans(cg_mean_exp_mat)
  # group wise sds
  cg_mean_cgw_sds <- apply(cg_mean_exp_mat, 1, sd)
  # group wise z-score
  cg_mean_cgw_zscore_mat <- sweep(cg_mean_exp_mat, 1, cg_mean_cgw_means,
                                  FUN = '-')
  cg_mean_cgw_zscore_mat <- sweep(cg_mean_cgw_zscore_mat, 1, cg_mean_cgw_sds,
                                  FUN = '/')

  cg_mean_cgw_zscore_df <- data.frame(
    t(cg_mean_cgw_zscore_mat), check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same
  stopifnot(identical(rownames(cg_mean_cgw_zscore_df), check_gids))
  # assert cancer groups are the same
  stopifnot(identical(sort(colnames(cg_mean_cgw_zscore_df)), check_groups))
  res_list$zscore_df <- cg_mean_cgw_zscore_df


  # compute ranks
  print('Compute quantiles...')
  # cg_mean_exp_df is (n_groups, n_genes)
  cg_mean_exp_mat <- as.matrix(cg_mean_exp_df)
  # assert gene symbols are the same
  stopifnot(identical(colnames(cg_mean_exp_mat), check_gids))
  # assert cancer groups are the same
  stopifnot(identical(sort(rownames(cg_mean_exp_mat)),
                      check_groups))
  # group wise ranks
  cg_mean_cgw_rank_mat <- apply(cg_mean_exp_mat, 1,
                                function(x) rank(x, ties.method='min'))
  # describes which quartile the genes are in
  cg_mean_cgw_d_mat <- cg_mean_cgw_rank_mat

  p0p25_idc <- cg_mean_cgw_rank_mat > (nrow(cg_mean_cgw_rank_mat) * 0.75)
  cg_mean_cgw_d_mat[p0p25_idc] <- 'Highest expressed 25%'

  p25p50_idc <- cg_mean_cgw_rank_mat > (nrow(cg_mean_cgw_rank_mat) * 0.5) &
    cg_mean_cgw_rank_mat <= (nrow(cg_mean_cgw_rank_mat) * 0.75)
  # paste0 to make the line shorter
  cg_mean_cgw_d_mat[p25p50_idc] <- paste0(
    'Expression between upper quartile and median')

  p50p75_idc <- cg_mean_cgw_rank_mat > (nrow(cg_mean_cgw_rank_mat) * 0.25) &
    cg_mean_cgw_rank_mat <= (nrow(cg_mean_cgw_rank_mat) * 0.5)
  cg_mean_cgw_d_mat[p50p75_idc] <- paste0(
    'Expression between median and lower quartile')

  p75p100_idc <- cg_mean_cgw_rank_mat <= (nrow(cg_mean_cgw_rank_mat) * 0.25)
  cg_mean_cgw_d_mat[p75p100_idc] <- 'Lowest expressed 25%'
  stopifnot(identical(
    sort(unique(as.vector(cg_mean_cgw_d_mat))),
    c("Expression between median and lower quartile",
      "Expression between upper quartile and median", 
      "Highest expressed 25%", "Lowest expressed 25%")
  ))

  cg_mean_cgw_quant_out_df <- data.frame(
    cg_mean_cgw_d_mat, check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same
  stopifnot(identical(rownames(cg_mean_cgw_quant_out_df), check_gids))
  # assert cancer groups are the same
  stopifnot(identical(sort(colnames(cg_mean_cgw_quant_out_df)),
                      check_groups))
  res_list$quant_df <- cg_mean_cgw_quant_out_df


  return(res_list)
}

# load tpm matrix and sample metadata -------------------------------------
# NOTE: TPM matrices will be replaced by batch effect removed TPM matrix
# in the future
# gtex target tcga gmkf pbta
tpm_df <- readRDS('../../data/gene-expression-rsem-tpm-collapsed.rds')

suppressWarnings(
  htl_df <- readr::read_tsv(
    '../../data/histologies.tsv',
    col_types = cols(
      Kids_First_Biospecimen_ID = col_character(),
      sample_id = col_character(),
      aliquot_id = col_character(),
      Kids_First_Participant_ID = col_character(),
      experimental_strategy = col_character(),
      sample_type = col_character(),
      composition = col_character(),
      tumor_descriptor = col_character(),
      primary_site = col_character(),
      reported_gender = col_character(),
      race = col_character(),
      ethnicity = col_character(),
      diagnosis_type = col_character(),
      diagnosis_category = col_character(),
      age_at_diagnosis_days = col_character(),
      pathology_diagnosis = col_character(),
      RNA_library = col_character(),
      EFS_days = col_double(),
      OS_days = col_double(),
      OS_status = col_character(),
      PFS_days = col_character(),
      cohort = col_character(),
      age_last_update_days = col_double(),
      seq_center = col_character(),
      parent_aliquot_id = col_character(),
      previous_parent_aliquot_id = col_character(),
      cancer_predispositions = col_character(),
      previous_cancer_predispositions = col_character(),
      pathology_free_text_diagnosis = col_character(),
      cohort_participant_id = col_character(),
      germline_sex_estimate = col_character(),
      extent_of_tumor_resection = col_character(),
      normal_fraction = col_double(),
      tumor_fraction = col_double(),
      tumor_ploidy = col_double(),
      CNS_region = col_character(),
      molecular_subtype = col_character(),
      integrated_diagnosis = col_character(),
      Notes = col_character(),
      harmonized_diagnosis = col_character(),
      broad_histology = col_character(),
      short_histology = col_character(),
      cancer_group = col_character(),
      gtex_group = col_character(),
      gtex_subgroup = col_character()
    )
  )
)
# assert read count matrix column names match metadata sample IDs
stopifnot(all(colnames(tpm_df) %in% htl_df$Kids_First_Biospecimen_ID))
# assert read count matrix column naems are unique
stopifnot(identical(
  ncol(tpm_df),
  length(unique(colnames(tpm_df)))
))
# assert metadata sample IDs are unique
stopifnot(identical(
  nrow(htl_df),
  length(unique(htl_df$Kids_First_Biospecimen_ID))
))

# subset histology data frame to have the same samples as the TPM matrix
rna_htl_df <- htl_df[htl_df$Kids_First_Biospecimen_ID %in% colnames(tpm_df), ]
stopifnot(identical(
  sort(rna_htl_df$Kids_First_Biospecimen_ID),
  sort(colnames(tpm_df))
))


# summary statistics of all cohorts --------------------------------------------
print('Compute TPM summary statistics for all cohorts...')
# select cancer samples
ac_rna_htl_df <- rna_htl_df[!is.na(rna_htl_df$cancer_group), ]
# select cancer_groups with >= 5 samples
cg_cnts <- table(ac_rna_htl_df$cancer_group)
nge5_cg_cnts <- cg_cnts[cg_cnts >= 5]
nge5_cgs <- names(nge5_cg_cnts)

c_rna_htl_df <- ac_rna_htl_df[ac_rna_htl_df$cancer_group %in% nge5_cgs, ]
c_tpm_df <- tpm_df[, c_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(c_tpm_df), rownames(tpm_df)))
# assert c_tpm_df rows match c_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  c_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(c_tpm_df)
))

# summary stats data frames
ot_tpm_ss_dfs <- get_expression_summary_stats(
  c_tpm_df, c_rna_htl_df$cancer_group)

write_tsv(rownames_to_column(ot_tpm_ss_dfs$mean_df, 'gene'),
          'results/cancer_group_all_cohort_mean_tpm.tsv')
write_tsv(rownames_to_column(ot_tpm_ss_dfs$sd_df, 'gene'),
          'results/cancer_group_all_cohort_standard_deviation_tpm.tsv')
write_tsv(
  rownames_to_column(ot_tpm_ss_dfs$zscore_df, 'gene'),
  'results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_z_scores.tsv')
write_tsv(
  rownames_to_column(ot_tpm_ss_dfs$quant_df, 'gene'),
  'results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_quantiles.tsv')



# summary statistics of each cohort ---------------------------------------
print('Compute TPM summary statistics for each cohort...')
# select cancer samples
acc_rna_htl_df <- rna_htl_df[!(is.na(rna_htl_df$cancer_group) |
                                is.na(rna_htl_df$cohort)), ]
acc_rna_htl_df$cancer_group_cohort <- paste(
  acc_rna_htl_df$cancer_group, acc_rna_htl_df$cohort, sep = '___')
# select cancer_group_cohorts with >= 5 samples
cgc_cnts <- table(acc_rna_htl_df$cancer_group_cohort)
nge5_cgc_cnts <- cgc_cnts[cgc_cnts >= 5]
nge5_cgcs <- names(nge5_cgc_cnts)

cc_rna_htl_df <- acc_rna_htl_df[
  acc_rna_htl_df$cancer_group_cohort %in% nge5_cgcs, ]
cc_tpm_df <- tpm_df[, cc_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(cc_tpm_df), rownames(tpm_df)))
# assert c_tpm_df rows match c_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  cc_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(cc_tpm_df)
))

# summary stats data frames
cc_tpm_ss_dfs <- get_expression_summary_stats(
  cc_tpm_df, cc_rna_htl_df$cancer_group_cohort)

write_tsv(rownames_to_column(cc_tpm_ss_dfs$mean_df, 'gene'),
          'results/cancer_group_individual_cohort_mean_tpm.tsv')
write_tsv(rownames_to_column(cc_tpm_ss_dfs$sd_df, 'gene'),
          'results/cancer_group_individual_cohort_standard_deviation_tpm.tsv')
write_tsv(
  rownames_to_column(cc_tpm_ss_dfs$zscore_df, 'gene'),
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_cancer_group_wise_mean_tpm_z_scores.tsv')))
write_tsv(
  rownames_to_column(cc_tpm_ss_dfs$quant_df, 'gene'),
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_cancer_group_wise_mean_tpm_quantiles.tsv')))
