suppressPackageStartupMessages(library(tidyverse))


# Generate sample metadata data frame.
# Args:
#   - htl_df: histology data frame. Columns must contain group_var and
#     Kids_First_Biospecimen_ID.
#   - group_var: character of the grouping variable, which is a column in the
#     htl_df.
#
# Returns a tibble of the number of samples and sample IDs of each group.
get_sample_meta_df <- function(htl_df, group_var) {
  mdf <- htl_df %>%
    group_by_at(group_var) %>%
    summarise(
      n_samples = n(),
      Kids_First_Biospecimen_IDs = paste(Kids_First_Biospecimen_ID,
                                         collapse=',')) %>%
    arrange_at(group_var)

  return(mdf)
}


# Get summary statistics dataframe for output.
#
# Args:
# - ss_df: (n_genes, n_groups) summary statistics data frame. Rownames must be
#   gene symbols.
# - gsb_gids_df: gene symbol and ENSG ID data frame. Two columns must be
#   gene_symbol and gene_id.
#
# Returns a (n_genes, n_groups) summary statistics tibble with first two
# columns as gene_id and gene_symbol.
get_output_ss_df <- function(ss_df, gsb_gid_df) {
  gsb_gids_conv_df <- gsb_gid_df
  # assert all symbols are unique
  stopifnot(identical(length(unique(gsb_gids_conv_df$gene_symbol)),
                      nrow(gsb_gids_conv_df)))
  rownames(gsb_gids_conv_df) <- gsb_gids_conv_df$gene_symbol
  # assert rownames are not changed automatically
  stopifnot(identical(rownames(gsb_gids_conv_df),
                      gsb_gids_conv_df$gene_symbol))

  stopifnot(all(rownames(ss_df) %in% rownames(gsb_gids_conv_df)))
  ss_gsb_gid_df <- gsb_gids_conv_df[rownames(ss_df), ]

  # assert rownames(ss_df) are the same as rownames(ss_gsb_gid_df)
  stopifnot(identical(rownames(ss_df), rownames(ss_gsb_gid_df)))
  ss_out_df <- cbind(ss_gsb_gid_df, ss_df)

  stopifnot(identical(rownames(ss_out_df), ss_out_df$gene_symbol))
  ss_out_df <- remove_rownames(ss_out_df)
  return(ss_out_df)
}


# Generate means, standard deviations, z-scores, and ranks within each group.
#
# Args:
# - exp_df: (n_genes, n_samples) expression level numeric data frame
# - groups: character vector of length n_samples, which is used for grouping
#   the samples.
#
# Returns a list of (n_genes, n_groups) summary statistics tables.
get_expression_summary_stats <- function(exp_df, groups) {
  # unique groups to check that the computing steps do not modify the groups.
  check_groups <- sort(unique(groups))
  # gene symbols to check that the computing steps do not modify
  # the rownames of the exp_df.
  check_gids <- rownames(exp_df)

  # set check.names = FALSE and check.rows = FALSE to avoid R from
  # changing the rownames or colnames implicitly
  c_exp_df <- data.frame(t(exp_df), check.names = FALSE, check.rows = FALSE)
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
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_mean_exp_out_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(cg_mean_exp_out_df)), check_groups))
  res_list$mean_df <- cg_mean_exp_out_df


  # compute group standard deviations
  print('Compute standard deviations...')
  cg_sd_exp_df <- c_exp_df %>%
    group_by(sample_group) %>%
    summarise_all(sd) %>%
    column_to_rownames('sample_group')

  cg_sd_exp_out_df <- data.frame(
    t(cg_sd_exp_df), check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_sd_exp_out_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(cg_sd_exp_out_df)), check_groups))
  res_list$sd_df <- cg_sd_exp_out_df

  # input is a numeric matrix
  # procedure adapted from @kgaonkar6's code at
  # <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/blob/
  #     0a85a711709b5adc1e56a26a397d238cb3ebbb58/analyses/
  #     fusion_filtering/03-Calc-zscore-annotate.R#L115-L121>
  row_wise_zscores <- function(num_mat) {
    row_means <- rowMeans(num_mat)
    row_sds <- apply(num_mat, 1, sd)
    # row-wise z-score
    row_wise_zscore_mat <- sweep(num_mat, 1, row_means, FUN = '-')
    row_wise_zscore_mat <- sweep(row_wise_zscore_mat, 1, row_sds, FUN = '/')

    return(row_wise_zscore_mat)
  }

  # compute z-scores
  print('Compute group-wise z-scores...')
  # cg_mean_exp_df is (n_groups, n_genes)
  cg_mean_exp_mat <- as.matrix(cg_mean_exp_df)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(colnames(cg_mean_exp_mat), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(rownames(cg_mean_exp_mat)), check_groups))
  # group wise z-scores
  cg_mean_cgw_zscore_mat <- row_wise_zscores(cg_mean_exp_mat)

  cg_mean_cgw_zscore_df <- data.frame(
    t(cg_mean_cgw_zscore_mat), check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_mean_cgw_zscore_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(cg_mean_cgw_zscore_df)), check_groups))
  res_list$group_wise_zscore_df <- cg_mean_cgw_zscore_df


  # compute z-scores
  print('Compute gene-wise z-scores...')
  # cg_mean_exp_df is (n_groups, n_genes)
  # so gr_cg_mean_exp_mat is (n_genes, n_groups)
  # gr = gene rows
  gr_cg_mean_exp_mat <- t(cg_mean_exp_df)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(gr_cg_mean_exp_mat), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(gr_cg_mean_exp_mat)), check_groups))
  cg_mean_gene_wise_zscore_mat <- row_wise_zscores(gr_cg_mean_exp_mat)

  cg_mean_gene_wise_zscore_df <- data.frame(
    cg_mean_gene_wise_zscore_mat, check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_mean_gene_wise_zscore_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(
    sort(colnames(cg_mean_gene_wise_zscore_df)), check_groups))
  res_list$gene_wise_zscore_df <- cg_mean_gene_wise_zscore_df


  # compute ranks
  print('Compute quantiles...')
  # cg_mean_exp_df is (n_groups, n_genes)
  cg_mean_exp_mat <- as.matrix(cg_mean_exp_df)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(colnames(cg_mean_exp_mat), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(rownames(cg_mean_exp_mat)), check_groups))
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
  # assert all entries have description values
  stopifnot(identical(
    sort(unique(as.vector(cg_mean_cgw_d_mat))),
    c("Expression between median and lower quartile",
      "Expression between upper quartile and median", 
      "Highest expressed 25%", "Lowest expressed 25%")
  ))

  cg_mean_cgw_quant_out_df <- data.frame(
    cg_mean_cgw_d_mat, check.names = FALSE, check.rows = FALSE)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_mean_cgw_quant_out_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(cg_mean_cgw_quant_out_df)), check_groups))
  res_list$quant_df <- cg_mean_cgw_quant_out_df


  return(res_list)
}

# load tpm matrix and sample metadata -------------------------------------
# NOTE: TPM matrices may be replaced by batch effect removed TPM matrix
# in the future
message('Read data...')

tpm_df <- readRDS('../../data/gene-expression-rsem-tpm-collapsed.rds')

htl_df <- readr::read_tsv(
  '../../data/histologies.tsv',
  col_types = cols(
    .default = col_guess(),
    gtex_group = col_character(),
    gtex_subgroup = col_character(),
    EFS_days = col_number()
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

# independent sample table
suppressMessages(
  rna_idps_df <- read_tsv(
    file.path('..', 'independent-samples', 'results',
              'independent-specimens.rnaseq.primary.tsv'))
)
rna_idps_kfbids <- rna_idps_df$Kids_First_Biospecimen_ID
# subset histology data frame to have the same samples as the TPM matrix
rna_htl_df <- htl_df[
  htl_df$Kids_First_Biospecimen_ID %in% colnames(tpm_df) &
    htl_df$Kids_First_Biospecimen_ID %in% rna_idps_kfbids, ]
stopifnot(all(rna_htl_df$Kids_First_Biospecimen_ID %in% colnames(tpm_df)))

# gene symbol to ENSG id table
gid_gsb_tbl <- read_tsv('../../data/ensg-hugo-rmtl-v1-mapping.tsv',
                        col_types = cols(.default = col_guess()),
                        guess_max = 10000) %>%
  rename(gene_id = ensg_id) %>%
  select(gene_id, gene_symbol) %>%
  distinct() 

# Collapse gid_gsb_tbl by gene_symbol.
# If one gene_symbol is mapped to multiple gene_ids (ENSG IDs), the value of
# the gene_id column in the collapsed gsb_gids_tbl is a comma separated list of
# all Ensembl ENSG IDs.
gsb_gids_tbl <- gid_gsb_tbl %>%
  group_by(gene_symbol) %>%
  summarise(gene_id = paste(gene_id, collapse = ','))

stopifnot(identical(length(unique(gsb_gids_tbl$gene_symbol)),
                    nrow(gsb_gids_tbl)))
gsb_gids_df <- data.frame(gsb_gids_tbl, stringsAsFactors = FALSE)


# summary statistics of all cohorts --------------------------------------------
print('Compute TPM summary statistics for all cohorts...')
# select cancer samples
ac_rna_htl_df <- rna_htl_df[!is.na(rna_htl_df$cancer_group), ]
# select cancer_groups with >= 5 samples
# adapted from https://stackoverflow.com/a/40091131/4638182
c_rna_htl_df <- ac_rna_htl_df %>%
  group_by(cancer_group) %>%
  filter(n() >= 5) %>%
  ungroup()

c_sample_meta_df <- get_sample_meta_df(c_rna_htl_df, 'cancer_group')
write_tsv(c_sample_meta_df,
          'results/cancer_group_all_cohort_sample_metadata.tsv')

c_tpm_df <- tpm_df[, c_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(c_tpm_df), rownames(tpm_df)))
# assert c_tpm_df columns match c_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  c_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(c_tpm_df)
))

# summary stats data frames
ot_tpm_ss_dfs <- get_expression_summary_stats(
  c_tpm_df, c_rna_htl_df$cancer_group)

write_tsv(get_output_ss_df(ot_tpm_ss_dfs$mean_df, gsb_gids_df),
          'results/cancer_group_all_cohort_mean_tpm.tsv')
write_tsv(get_output_ss_df(ot_tpm_ss_dfs$sd_df, gsb_gids_df),
          'results/cancer_group_all_cohort_standard_deviation_tpm.tsv')
write_tsv(
  get_output_ss_df(ot_tpm_ss_dfs$group_wise_zscore_df, gsb_gids_df),
  'results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv')
write_tsv(
  get_output_ss_df(ot_tpm_ss_dfs$quant_df, gsb_gids_df),
  'results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv')
write_tsv(
  get_output_ss_df(ot_tpm_ss_dfs$gene_wise_zscore_df, gsb_gids_df),
  'results/cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv')



# summary statistics of each cohort ---------------------------------------
print('Compute TPM summary statistics for each cohort...')
# select cancer samples
acc_rna_htl_df <- rna_htl_df[!(is.na(rna_htl_df$cancer_group) |
                                is.na(rna_htl_df$cohort)), ]
acc_rna_htl_df$cancer_group_cohort <- paste(
  acc_rna_htl_df$cancer_group, acc_rna_htl_df$cohort, sep = '___')
# select cancer_group_cohorts with >= 5 samples
# adapted from https://stackoverflow.com/a/40091131/4638182
cc_rna_htl_df <- acc_rna_htl_df %>%
  group_by(cancer_group_cohort) %>%
  filter(n() >= 5) %>%
  ungroup()

cc_sample_meta_df <- get_sample_meta_df(cc_rna_htl_df, 'cancer_group_cohort')
write_tsv(cc_sample_meta_df,
          'results/cancer_group_individual_cohort_sample_metadata.tsv')


cc_tpm_df <- tpm_df[, cc_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(cc_tpm_df), rownames(tpm_df)))
# assert cc_tpm_df columns match cc_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  cc_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(cc_tpm_df)
))

# summary stats data frames
cc_tpm_ss_dfs <- get_expression_summary_stats(
  cc_tpm_df, cc_rna_htl_df$cancer_group_cohort)

write_tsv(get_output_ss_df(cc_tpm_ss_dfs$mean_df, gsb_gids_df),
          'results/cancer_group_individual_cohort_mean_tpm.tsv')
write_tsv(get_output_ss_df(cc_tpm_ss_dfs$sd_df, gsb_gids_df),
          'results/cancer_group_individual_cohort_standard_deviation_tpm.tsv')
write_tsv(
  get_output_ss_df(cc_tpm_ss_dfs$group_wise_zscore_df, gsb_gids_df),
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_cancer_group_wise_z_scores.tsv')))
write_tsv(
  get_output_ss_df(cc_tpm_ss_dfs$quant_df, gsb_gids_df),
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_cancer_group_wise_quantiles.tsv')))
write_tsv(
  get_output_ss_df(cc_tpm_ss_dfs$gene_wise_zscore_df, gsb_gids_df),
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_gene_wise_z_scores.tsv')))



# Generate combined long table -------------------------------------------------
# - Combine cancer_group and cancer_group_cohort TPM mean/SD/zscore/quantile
#   into one table. 
# - One type of z-score per table.
# - Row-wise: Ensembl gene ID, gene symbol, cancer_group_cohort, TPM mean,
#   TPM SD, TPM mean zscore, TPM mean quantile.

# Convert a list of wide summary stat dfs into a single long df
#
# Args:
# - sum_stat_df_list: a list of (n_genes, n_groups) summary statistics data
#   frame. Rownames must be gene symbols. The names of sum_stat_df_list must in
#   c("mean_df", "sd_df", "group_wise_zscore_df", "gene_wise_zscore_df", 
#   "quant_df"). Must contain "mean_df" table.
# - gsb_gids_df: gene symbol and ENSG ID data frame. Two columns must be
#   gene_symbol and gene_id. Passed to get_output_ss_df.
# - key_colname: column name in the converted long df that represents the
#   columns of the wide data frame.
# - sample_meta_df: sample information of each column in the wide tables,
#   generated by get_sample_meta_df
#
# Returns a long-format tibble with each row has Ensembl gene ID, gene symbol,
# cancer_group_cohort, TPM mean, TPM SD, TPM mean zscore, TPM mean quantile.
#
# Note: If one gene symbol matches to multiple Ensembl gene IDs, each Ensembl
# gene ID will become one row in the result table, so gene symbols and summary
# statistics have duplicated rows.
sum_stat_df_list_wide_to_long <- function(sum_stat_df_list, gsb_gids_df,
                                          key_colname, sample_meta_df) {
  sum_stats_df_names_vec <- names(sum_stat_df_list)
  names(sum_stats_df_names_vec) <- sum_stats_df_names_vec
  
  long_tbls <- lapply(sum_stats_df_names_vec, function(x) {
    x_df <- sum_stat_df_list[[x]]
    if (x == 'mean_df') {
      val_colname <- 'tpm_mean'
    } else if (x == 'sd_df') {
      val_colname <- 'tpm_sd'
    } else if (x == 'group_wise_zscore_df') {
      val_colname <- 'tpm_mean_cancer_group_wise_zscore'
    } else if (x == 'gene_wise_zscore_df') {
      val_colname <- 'tpm_mean_gene_wise_zscore'
    } else if (x == 'quant_df') {
      val_colname <- 'tpm_mean_cancer_group_wise_quantiles'
    } else {
      stop(paste('unknown df', x))
    }

    # the !! syntax is found at https://stackoverflow.com/a/54013082/4638182
    # pivot_longer is not available for the docker image with dplyr 0.8.x
    x_tbl <- get_output_ss_df(x_df, gsb_gids_df) %>%
      separate_rows(gene_id, sep = ',') %>%
      gather(key = !!key_colname, value = !!val_colname, -gene_symbol, -gene_id)
    return(x_tbl)
  })

  m_long_tbl <- reduce(long_tbls, function(x, y) {
    return(full_join(x, y, by = c('gene_symbol', 'gene_id', key_colname)))
  })
  m_ann_long_tbl <- m_long_tbl %>%
    left_join(sample_meta_df, by = key_colname)
  stopifnot(identical(
    sum(is.na(m_ann_long_tbl[, c('gene_symbol', 'gene_id', key_colname,
                                 'n_samples', 'Kids_First_Biospecimen_IDs')])),
    as.integer(0)))
  stopifnot(identical(nrow(m_ann_long_tbl), nrow(m_long_tbl)))
  stopifnot(identical(sort(unique(m_ann_long_tbl$gene_symbol)),
                      sort(unique(rownames(sum_stat_df_list$mean_df)))))
  stopifnot(identical(sort(unique(m_ann_long_tbl[, key_colname, drop = TRUE])),
                      sort(colnames(sum_stat_df_list$mean_df))))
  return(m_ann_long_tbl)
}

ot_tpm_ss_long_tbl <- sum_stat_df_list_wide_to_long(
  ot_tpm_ss_dfs, gsb_gids_df, 'cancer_group', c_sample_meta_df) %>%
  mutate(cancer_group_cohort = paste0(cancer_group, '___AllCohorts')) %>%
  select(-cancer_group) %>%
  select(gene_symbol, gene_id, cancer_group_cohort, tpm_mean, 
         tpm_sd, tpm_mean_cancer_group_wise_zscore, tpm_mean_gene_wise_zscore, 
         tpm_mean_cancer_group_wise_quantiles, n_samples,
         Kids_First_Biospecimen_IDs)

cc_tpm_ss_long_tbl <- sum_stat_df_list_wide_to_long(
  cc_tpm_ss_dfs, gsb_gids_df, 'cancer_group_cohort', cc_sample_meta_df)

stopifnot(identical(colnames(ot_tpm_ss_long_tbl),
                    colnames(cc_tpm_ss_long_tbl)))

m_tpm_ss_long_tbl <- bind_rows(ot_tpm_ss_long_tbl, cc_tpm_ss_long_tbl)
stopifnot(identical(colnames(m_tpm_ss_long_tbl),
                    colnames(cc_tpm_ss_long_tbl)))
m_tpm_ss_long_tbl <- m_tpm_ss_long_tbl %>%
  separate(col = cancer_group_cohort, into = c('cancer_group', 'cohort'),
           sep = '___')
stopifnot(identical(nrow(distinct(m_tpm_ss_long_tbl)), 
                    nrow(m_tpm_ss_long_tbl)))
stopifnot(identical(
  sum(is.na(m_tpm_ss_long_tbl[, c('gene_symbol', 'gene_id', 'cancer_group',
                                  'cohort', 'n_samples',
                                  'Kids_First_Biospecimen_IDs')])),
  as.integer(0)))

write_tsv(select(m_tpm_ss_long_tbl,
                 -tpm_mean_gene_wise_zscore,
                 -Kids_First_Biospecimen_IDs),
          'results/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv')

write_tsv(select(m_tpm_ss_long_tbl,
                 -tpm_mean_cancer_group_wise_zscore,
                 -Kids_First_Biospecimen_IDs),
          'results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv')
