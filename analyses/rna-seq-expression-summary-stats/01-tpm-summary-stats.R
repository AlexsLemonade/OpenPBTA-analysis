suppressPackageStartupMessages(library(tidyverse))



# Function definitions ---------------------------------------------------------
# Favor function definitions in one file with data processing code over `source`
# other scripts. `source` is not explicit on the function name of the imported
# function. Even though the script can be named the same as the funtion name,
# it is unclear whether other functions are also imported.


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
    cg_mean_cgw_d_mat, check.names = FALSE, check.rows = FALSE,
    stringsAsFactors = FALSE)
  # assert gene symbols are the same as input data frame
  stopifnot(identical(rownames(cg_mean_cgw_quant_out_df), check_gids))
  # assert groups are the same as input data frame
  stopifnot(identical(sort(colnames(cg_mean_cgw_quant_out_df)), check_groups))
  res_list$quant_df <- cg_mean_cgw_quant_out_df


  return(res_list)
}


# Get summary statistics dataframe for output.
#
# Args:
# - ss_df: (n_genes, n_groups) summary statistics data frame. Rownames must be
#   gene symbols.
# - gsb_gid_df: gene symbol and ENSG ID data frame. Two columns must be
#   gene_symbol and gene_id. ENSG ID column can be a comma separated list of
#   ENSG IDs for gene symbols that are mapped to multiple ENSG IDs.
#
# Returns a (n_genes, n_groups + 2) summary statistics tibble with first two
# columns as gene_id and gene_symbol.
get_output_ss_df <- function(ss_df, gsb_gid_df) {
  # Add gene_id and gene_symbol to each row
  # gene_symbol to gene_ids table for annotation
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

  ss_out_df <- as_tibble(remove_rownames(ss_out_df))
  # assert colnames are not changed
  stopifnot(identical(colnames(ss_df),
                      colnames(select(ss_out_df, -gene_id, -gene_symbol))))

  # Replace NA/NaNs with ''
  # assert no NA/NaNs in gene_id or gene_symbol
  stopifnot(identical(
    as.integer(0),
    sum(is.na(select(ss_out_df, gene_id, gene_symbol)))
  ))
  rm_na_ss_out_df <- mutate_all(
    ss_out_df, function(x) replace_na(x, replace = ''))
  stopifnot(identical(as.integer(0), sum(is.na(rm_na_ss_out_df))))
  return(rm_na_ss_out_df)
}


# Generate means, standard deviations, group/gene-wise z-scores, and ranks
# within each group for output
#
# Args:
# - exp_df: (n_genes, n_samples) expression level numeric data frame
# - groups: character vector of length n_samples, which is used for grouping the
#   samples.
# - gsb_gid_df: gene symbol and ENSG ID data frame. Two columns must be
#   gene_symbol and gene_id. ENSG ID column can be a comma separated list of
#   ENSG IDs for gene symbols that are mapped to multiple ENSG IDs.
#
# Returns a list of (n_genes, n_groups) summary statistics tibbles with first
# two columns as gene_id and gene_symbol.
get_expression_summary_stats_out_dfs <- function(exp_df, groups, gsb_gid_df) {
  exp_ss_dfs <- get_expression_summary_stats(exp_df, groups)
  exp_ss_out_dfs <- lapply(exp_ss_dfs, function(x) {
    return(get_output_ss_df(x, gsb_gid_df))
  })
  return(exp_ss_out_dfs)
}


# Convert a list of wide summary stat dfs into a single long df
#
# Args:
# - sum_stat_out_df_list: a list of (n_genes, n_groups + 2) (+2 are gene_id and
#   gene_symbol) summary statistics output data frame. Columns must contain
#   gene_symbol and gene_id. The names of sum_stat_out_df_list must in
#   c("mean_df", "sd_df", "group_wise_zscore_df", "gene_wise_zscore_df",
#   "quant_df"). Must contain "mean_df" table.
# - key_colname: column name in the converted long df that represents the
#   columns of the wide data frame.
# - sample_meta_df: sample information of each column in the wide tables. Must
#   contain key_colname to join by.
#
# Returns a long-format tibble with each row has Ensembl gene ID, gene symbol,
# cancer_group_cohort, TPM mean, TPM SD, TPM mean zscore, TPM mean quantile.
#
# Note: If one gene symbol matches to multiple Ensembl gene IDs, each Ensembl
# gene ID will become one row in the result table, so gene symbols and summary
# statistics have duplicated rows.
sum_stat_df_list_wide_to_long <- function(sum_stat_out_df_list, key_colname,
                                          sample_meta_df) {
  sum_stats_df_names_vec <- names(sum_stat_out_df_list)
  names(sum_stats_df_names_vec) <- sum_stats_df_names_vec

  long_tbls <- lapply(sum_stats_df_names_vec, function(x) {
    x_df <- sum_stat_out_df_list[[x]]
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
    x_tbl <- x_df %>%
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
    sum(is.na(m_ann_long_tbl[, c('gene_symbol', 'gene_id',
                                 colnames(sample_meta_df))])),
    as.integer(0)))
  stopifnot(identical(nrow(m_ann_long_tbl), nrow(m_long_tbl)))
  stopifnot(identical(sort(unique(m_ann_long_tbl$gene_symbol)),
                      sort(unique(sum_stat_out_df_list$mean_df$gene_symbol))))
  stopifnot(identical(sort(unique(m_ann_long_tbl[, key_colname, drop = TRUE])),
                      sort(colnames(sum_stat_out_df_list$mean_df %>%
                                      select(-gene_symbol, -gene_id)))))
  return(m_ann_long_tbl)
}


# Formats a character vector of cohorts to
# cohort1_n1_samples[&cohort2_n2_samples&cohort3_n3_samples]
#
# Args:
# - cohort_vec: a character vector of cohorts
#
# Returns a string in the format of
# cohort1_n1_samples[&cohort2_n2_samples&cohort3_n3_samples]
format_cohort_sample_counts <- function(cohort_vec) {
  stopifnot(identical(as.integer(0), sum(is.na(cohort_vec))))
  stopifnot(is.character(cohort_vec))

  cohort_df <- tibble(cohort = cohort_vec) %>%
    group_by(cohort) %>%
    summarise(cohort_n_samples = paste(unique(cohort),
                                       length(cohort),
                                       'samples', sep = '_')) %>%
    arrange(cohort)

  return(paste(cohort_df$cohort_n_samples, collapse = '&'))
}
# # test case
# format_cohort_sample_counts(c('b', 'a', 'a', 'c', 'c'))
# format_cohort_sample_counts(c('a', 'a', 'b'))
# format_cohort_sample_counts(c('a', 'a'))
# format_cohort_sample_counts(c('a'))
# format_cohort_sample_counts(character(0))
# # following cases should fail
# format_cohort_sample_counts(c('a', NA))
# format_cohort_sample_counts(c(NA))
# format_cohort_sample_counts(c(NaN))



# Create output dir ------------------------------------------------------------
tables_dir <- 'results'

if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)
}



# Read data --------------------------------------------------------------------
# NOTE: TPM matrices may be replaced by batch effect removed TPM matrix
# in the future
message('Read data...')

tpm_df <- readRDS('../../data/gene-expression-rsem-tpm-collapsed.rds')
# assert read count matrix column naems are unique
stopifnot(identical(
  ncol(tpm_df),
  length(unique(colnames(tpm_df)))
))

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
# assert metadata sample IDs have no NA
stopifnot(identical(
  as.integer(0),
  sum(is.na(htl_df$Kids_First_Biospecimen_ID))
))
# assert metadata cohort have no NA
stopifnot(identical(as.integer(0), sum(is.na(htl_df$cohort))))
# assert metadata sample IDs are unique
stopifnot(identical(
  nrow(htl_df),
  length(unique(htl_df$Kids_First_Biospecimen_ID))
))


# independent sample table
rna_idps_df <- read_tsv(
  file.path('..', 'independent-samples', 'results',
            'independent-specimens.rnaseq.primary.eachcohort.tsv'),
  col_types = cols(.default = col_guess()))

rna_idps_kfbids <- rna_idps_df$Kids_First_Biospecimen_ID
# assert no NA in rna_idps_kfbids
stopifnot(identical(as.integer(0), sum(is.na(rna_idps_kfbids))))
# assert all rna_idps_kfbids are unique
stopifnot(identical(length(rna_idps_kfbids), length(unique(rna_idps_kfbids))))


# subset histology data frame to have independent samples in tpm_df
rna_htl_df <- htl_df %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(tpm_df),
         Kids_First_Biospecimen_ID %in% rna_idps_kfbids)
stopifnot(all(rna_htl_df$Kids_First_Biospecimen_ID %in% colnames(tpm_df)))


# gene symbol to ENSG id table
gid_gsb_tbl <- read_tsv('../../data/ensg-hugo-rmtl-v1-mapping.tsv',
                        col_types = cols(.default = col_guess()),
                        guess_max = 10000) %>%
  rename(gene_id = ensg_id) %>%
  select(gene_id, gene_symbol) %>%
  distinct()
# assert there is no , in gene_id
stopifnot(identical(
  as.integer(0),
  sum(is.na(str_detect(gid_gsb_tbl$gene_id, ',')))
))
# assert all gene_id and gene_symbols are not NA
stopifnot(identical(sum(is.na(gid_gsb_tbl$gene_id)), as.integer(0)))
stopifnot(identical(sum(is.na(gid_gsb_tbl$gene_symbol)), as.integer(0)))


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


# Summary statistics of all cohorts --------------------------------------------
message('Compute TPM summary statistics for all cohorts...')
# select cancer samples
# nf = n_samples filtered
nf_all_cohorts_rna_htl_df <- rna_htl_df %>%
  filter(!is.na(cancer_group), !is.na(cohort)) %>%
  group_by(cancer_group) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  mutate(sample_group = paste0(cancer_group, '___all_cohorts'))
# The sample_group column needs to describe three things of a sample:
# - cancer_group
# - cohort(s)
# - The way all samples are grouped together, either by all-cohorts or
#   each-individual-cohort in this analysis. The gene-wise z-scores for the same
#   (cancer_group, cohort) will be different when computed using all-cohorts or
#   each-cohort grouping method, if other samples are not identically grouped
#   identically in all-cohort or each-cohort method.
#
# For easier implementation:
# - cohort(s) are implicit in all-cohorts grouping method
# - grouping method is implicit in each-cohort grouping method


nf_all_cohorts_sample_meta_out_df <- nf_all_cohorts_rna_htl_df %>%
  group_by(sample_group) %>%
  summarise(
    n_samples = n(),
    Kids_First_Biospecimen_IDs = paste(Kids_First_Biospecimen_ID, collapse=','),
    included_one_or_more_cohorts = format_cohort_sample_counts(cohort)) %>%
  select(sample_group, included_one_or_more_cohorts, n_samples,
         Kids_First_Biospecimen_IDs) %>%
  arrange(sample_group)

write_tsv(nf_all_cohorts_sample_meta_out_df,
          'results/cancer_group_all_cohort_sample_metadata.tsv')

# This table is used for annotating long tables in the last section.
nf_all_cohorts_sample_meta_df <- nf_all_cohorts_rna_htl_df %>%
  group_by(sample_group) %>%
  summarise(
    n_samples = n(),
    cancer_group = unique(cancer_group),
    cohort = 'all_cohorts') %>%
  arrange(sample_group)


nf_all_cohorts_tpm_df <- tpm_df[
  , nf_all_cohorts_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(nf_all_cohorts_tpm_df), rownames(tpm_df)))
# assert nf_all_cohorts_tpm_df columns match
# nf_all_cohorts_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  nf_all_cohorts_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(nf_all_cohorts_tpm_df)
))

# summary stats data frames
nf_all_cohorts_tpm_ss_out_dfs <- get_expression_summary_stats_out_dfs(
  nf_all_cohorts_tpm_df, nf_all_cohorts_rna_htl_df$sample_group,
  gsb_gids_df)

write_tsv(nf_all_cohorts_tpm_ss_out_dfs$mean_df,
          'results/cancer_group_all_cohort_mean_tpm.tsv')
write_tsv(nf_all_cohorts_tpm_ss_out_dfs$sd_df,
          'results/cancer_group_all_cohort_standard_deviation_tpm.tsv')
write_tsv(
  nf_all_cohorts_tpm_ss_out_dfs$group_wise_zscore_df,
  'results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv')
write_tsv(
  nf_all_cohorts_tpm_ss_out_dfs$quant_df,
  'results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv')
write_tsv(
  nf_all_cohorts_tpm_ss_out_dfs$gene_wise_zscore_df,
  'results/cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv')



# Summary statistics of each cohort --------------------------------------------
message('Compute TPM summary statistics for each individual cohort...')
# select cancer samples
# ind = individual
nf_ind_cohort_rna_htl_df <- rna_htl_df %>%
  filter(!is.na(cancer_group), !is.na(cohort)) %>%
  mutate(sample_group = paste(cancer_group, cohort, sep = '___')) %>%
  group_by(sample_group) %>%
  filter(n() >= 5) %>%
  ungroup()

nf_ind_cohort_sample_meta_out_df <- nf_ind_cohort_rna_htl_df %>%
  group_by(sample_group) %>%
  summarise(
    n_samples = n(),
    Kids_First_Biospecimen_IDs = paste(Kids_First_Biospecimen_ID,
                                       collapse=',')) %>%
  arrange(sample_group)

write_tsv(nf_ind_cohort_sample_meta_out_df,
          'results/cancer_group_individual_cohort_sample_metadata.tsv')

# This table is used for annotating long tables in the last section.
nf_ind_cohort_sample_meta_df <- nf_ind_cohort_rna_htl_df %>%
  group_by(sample_group) %>%
  summarise(
    n_samples = n(),
    cancer_group = unique(cancer_group),
    cohort = unique(cohort)) %>%
  arrange(sample_group)

nf_ind_cohort_tpm_df <- tpm_df[
  , nf_ind_cohort_rna_htl_df$Kids_First_Biospecimen_ID]
# assert gene symbols are the same
stopifnot(identical(rownames(nf_ind_cohort_tpm_df), rownames(tpm_df)))
# assert nf_ind_cohort_tpm_df columns match
# nf_ind_cohort_rna_htl_df$Kids_First_Biospecimen_ID
stopifnot(identical(
  nf_ind_cohort_rna_htl_df$Kids_First_Biospecimen_ID,
  colnames(nf_ind_cohort_tpm_df)
))

# summary stats data frames
nf_ind_cohort_tpm_ss_out_dfs <- get_expression_summary_stats_out_dfs(
  nf_ind_cohort_tpm_df, nf_ind_cohort_rna_htl_df$sample_group,
  gsb_gids_df)

write_tsv(nf_ind_cohort_tpm_ss_out_dfs$mean_df,
          'results/cancer_group_individual_cohort_mean_tpm.tsv')
write_tsv(nf_ind_cohort_tpm_ss_out_dfs$sd_df,
          'results/cancer_group_individual_cohort_standard_deviation_tpm.tsv')
write_tsv(
  nf_ind_cohort_tpm_ss_out_dfs$group_wise_zscore_df,
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_cancer_group_wise_z_scores.tsv')))
write_tsv(
  nf_ind_cohort_tpm_ss_out_dfs$quant_df,
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_cancer_group_wise_quantiles.tsv')))
write_tsv(
  nf_ind_cohort_tpm_ss_out_dfs$gene_wise_zscore_df,
  file.path('results',
            paste0('cancer_group_individual_cohort',
                   '_mean_tpm_gene_wise_z_scores.tsv')))



# Generate combined long table -------------------------------------------------
# - Combine all-cohorts and each-cohort TPM mean/SD/zscore/quantile
#   into one table.
# - One type of z-score per table.
# - Row-wise: Ensembl gene ID, gene symbol, cancer_group, cohort, TPM mean,
#   TPM SD, TPM mean zscore, TPM mean quantile.
message('Generate combined long table.')
nf_all_cohorts_tpm_ss_long_tbl <- sum_stat_df_list_wide_to_long(
  nf_all_cohorts_tpm_ss_out_dfs, 'sample_group',
  nf_all_cohorts_sample_meta_df)

nf_ind_cohort_tpm_ss_long_tbl <- sum_stat_df_list_wide_to_long(
  nf_ind_cohort_tpm_ss_out_dfs, 'sample_group',
  nf_ind_cohort_sample_meta_df)

stopifnot(identical(colnames(nf_all_cohorts_tpm_ss_long_tbl),
                    colnames(nf_ind_cohort_tpm_ss_long_tbl)))

m_tpm_ss_long_tbl <- bind_rows(nf_all_cohorts_tpm_ss_long_tbl,
                               nf_ind_cohort_tpm_ss_long_tbl)
stopifnot(identical(colnames(m_tpm_ss_long_tbl),
                    colnames(nf_ind_cohort_tpm_ss_long_tbl)))

m_tpm_ss_long_tbl <- m_tpm_ss_long_tbl %>%
  select(gene_symbol, gene_id, cancer_group, cohort, tpm_mean, tpm_sd,
         tpm_mean_cancer_group_wise_zscore, tpm_mean_gene_wise_zscore,
         tpm_mean_cancer_group_wise_quantiles, n_samples)

stopifnot(identical(nrow(distinct(m_tpm_ss_long_tbl)),
                    nrow(m_tpm_ss_long_tbl)))
stopifnot(identical(
  sum(is.na(m_tpm_ss_long_tbl[, c('gene_symbol', 'gene_id', 'cancer_group',
                                  'cohort', 'n_samples')])),
  as.integer(0)))



# Add EFO, MONDO, RMTL to long tables ------------------------------------------
ann_efo_mondo_cg_df <- efo_mondo_cg_df %>%
  rename(EFO = efo_code, MONDO = mondo_code)

m_tpm_ss_long_tbl <- m_tpm_ss_long_tbl %>%
  left_join(ann_efo_mondo_cg_df, by = 'cancer_group') %>%
  replace_na(list(EFO = '', MONDO = ''))
stopifnot(identical(sum(is.na(m_tpm_ss_long_tbl)), as.integer(0)))

# asert all rmtl NAs have version NAs, vice versa
stopifnot(identical(is.na(ensg_hugo_rmtl_df$rmtl),
                    is.na(ensg_hugo_rmtl_df$version)))
ann_ensg_hugo_rmtl_df <- ensg_hugo_rmtl_df %>%
  select(ensg_id, rmtl, version) %>%
  filter(!is.na(rmtl), !is.na(version)) %>%
  mutate(RMTL = paste0(rmtl, ' (', version, ')')) %>%
  select(ensg_id, RMTL) %>%
  rename(gene_id = ensg_id)

m_tpm_ss_long_tbl <- m_tpm_ss_long_tbl %>%
  left_join(ann_ensg_hugo_rmtl_df, by = 'gene_id') %>%
  replace_na(list(RMTL = ''))
stopifnot(identical(sum(is.na(m_tpm_ss_long_tbl)), as.integer(0)))

m_tpm_ss_long_tbl <- m_tpm_ss_long_tbl %>%
  select(gene_symbol, RMTL, gene_id,
         cancer_group, EFO, MONDO, n_samples, cohort,
         tpm_mean, tpm_sd,
         tpm_mean_cancer_group_wise_zscore, tpm_mean_gene_wise_zscore,
         tpm_mean_cancer_group_wise_quantiles)


# Output long tables -----------------------------------------------------------
# Output:
# - one table that only contains gene-wise zscores
# - one table that only contains group-wise zscores
group_wise_zscore_m_tpm_ss_long_tbl <- select(
  m_tpm_ss_long_tbl, -tpm_mean_gene_wise_zscore)

gene_wise_zscore_m_tpm_ss_long_tbl <- select(
  m_tpm_ss_long_tbl, -tpm_mean_cancer_group_wise_zscore)

write_tsv(
  group_wise_zscore_m_tpm_ss_long_tbl,
  file.path(tables_dir, 'long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv'))

write_tsv(
  gene_wise_zscore_m_tpm_ss_long_tbl,
  file.path(tables_dir, 'long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv'))

jsonlite::write_json(
  group_wise_zscore_m_tpm_ss_long_tbl,
  file.path(tables_dir, 'long_n_tpm_mean_sd_quantile_group_wise_zscore.json'))

jsonlite::write_json(
  gene_wise_zscore_m_tpm_ss_long_tbl,
  file.path(tables_dir, 'long_n_tpm_mean_sd_quantile_gene_wise_zscore.json'))


message('Done running 01-tpm-summary-stats.R.')
