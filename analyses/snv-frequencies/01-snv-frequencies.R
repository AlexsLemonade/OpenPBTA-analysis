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
  fxv <- unique(reduce(fxl, c, .init = character(0)))
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
# - maf_df: a MAF tibble. Must contain the following fields:
#   Kids_First_Biospecimen_ID,
#   Variant_ID,
#   Gene_symbol,
#   Gene_full_name,
#   dbSNP_ID,
#   VEP_impact,
#   SIFT_impact,
#   PolyPhen_impact,
#   Variant_classification,
#   Variant_type,
#   mRNA_RefSeq_ID,
#   Protein_RefSeq_ID,
#   Gene_Ensembl_ID,
#   Protein_Ensembl_ID,
#   Protein_change,
#   HotSpot.
# - overall_histology_df: the histology tibble that contains all samples. Must
#   contain the following fields: Kids_First_Biospecimen_ID,
#   Kids_First_Participant_ID, cancer_group, and cohort. This
#   is used for computing the following columns: Total_mutations,
#   Patients_in_dataset, Total_mutations_Over_Patients_in_dataset and
#   Frequency_in_overall_dataset.
# - primary_histology_df: the histology tibble that contains primary tumor
#   samples. Must contain the Kids_First_Biospecimen_ID field.
# - relapse_histology_df: the histology tibble that contains relapse tumor
#   samples.
#   Must contain the Kids_First_Biospecimen_ID field.
# - ss_cancer_group: a character value of the cancer group to compute for.
# - ss_cohorts: a vector of character values of the cohorts to compute for.
#
# Returns a MAF tibble with additional frequency columns.
get_cg_ch_mut_freq_tbl <- function(maf_df, overall_histology_df,
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


  ss_maf_df <- maf_df %>%
    filter(Kids_First_Biospecimen_ID %in% ss_htl_df$Kids_First_Biospecimen_ID)

  # gc(reset = TRUE) If one Variant_ID has more than one records in one
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
  sample_var_df <- ss_maf_df %>%
    select(Kids_First_Biospecimen_ID, Variant_ID) %>%
    distinct()

  patient_var_df <- sample_var_df %>%
    left_join(bpid_ss_htl_df, by = 'Kids_First_Biospecimen_ID') %>%
    group_by(Variant_ID) %>%
    summarise(Total_mutations = length(unique(Kids_First_Participant_ID))) %>%
    mutate(Patients_in_dataset = ss_n_patients) %>%
    mutate(Total_mutations_Over_Patients_in_dataset =
             paste(Total_mutations, Patients_in_dataset, sep = '/')) %>%
    mutate(Frequency_in_overall_dataset = num_to_pct_chr(
      Total_mutations / Patients_in_dataset))

  # td = tumor descriptor
  td_var_df <- sample_var_df %>%
    group_by(Variant_ID) %>%
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

  # If one Variant_ID has more than 1 values in the summarised fields, add
  # code to handle duplicates.
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
              mRNA_RefSeq_ID = unique(RefSeq),
              Protein_RefSeq_ID = unique(Protein_RefSeq_ID),
              Gene_Ensembl_ID = unique(Gene),
              Protein_Ensembl_ID = unique(ENSP),
              Protein_change = unique(HGVSp_Short),
              HotSpotAllele = unique(HotSpotAllele)) %>%
    left_join(patient_var_df, by = 'Variant_ID') %>%
    left_join(td_var_df, by = 'Variant_ID') %>%
    mutate(Disease = ss_cancer_group,
           Dataset = paste(ss_cohorts, collapse = '&'),
           HotSpot = if_else(HotSpotAllele == 1, true = 'Y', false = 'N')) %>%
    arrange(desc(as.numeric(Total_mutations))) %>%
    select(Gene_symbol, Dataset, Disease, Variant_ID, dbSNP_ID,
           VEP_impact, SIFT_impact, PolyPhen_impact,Variant_classification,
           Variant_type, Gene_full_name, Protein_RefSeq_ID,
           Gene_Ensembl_ID, Protein_Ensembl_ID, Protein_change,
           Total_mutations_Over_Patients_in_dataset,
           Frequency_in_overall_dataset,
           Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
           Frequency_in_primary_tumors,
           Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
           Frequency_in_relapse_tumors, HotSpot) %>%
    replace_na(list(
      Gene_symbol = '', Dataset = '', Disease = '',
      Variant_ID = '', dbSNP_ID = '', VEP_impact = '',
      SIFT_impact = '', PolyPhen_impact = '',
      Variant_classification = '', Variant_type = '',
      Gene_full_name = '', Protein_RefSeq_ID = '',
      Gene_Ensembl_ID = '', Protein_Ensembl_ID = '',
      Protein_change = '',
      Total_mutations_Over_Patients_in_dataset = '',
      Frequency_in_overall_dataset = '',
      Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset = '',
      Frequency_in_primary_tumors = '',
      Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset = '',
      Frequency_in_relapse_tumors = '', HotSpot = ''))

  return(output_var_df)
}



# Create output dir ------------------------------------------------------------
tables_dir <- 'results'

if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)
}



# Read SNV, independent sample list, and histology data ------------------------
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
  '../independent-samples/results/independent-specimens.wgs.primary.tsv',
  col_types = cols(
    .default = col_guess()))

# relapse independent samples
relapse_indp_sdf <- read_tsv(
  '../independent-samples/results/independent-specimens.wgs.relapse.tsv',
  col_types = cols(
    .default = col_guess()))



# Subset independent samples in histology table --------------------------------
message('Prepare data...')
# td = tumor descriptor
td_htl_dfs <- list(
  primary_htl_df = htl_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             primary_indp_sdf$Kids_First_Biospecimen_ID),

  relapse_htl_df = htl_df %>%
    filter(Kids_First_Biospecimen_ID %in%
             relapse_indp_sdf$Kids_First_Biospecimen_ID),

  overall_htl_df = htl_df %>%
    filter(sample_type == 'Tumor')
)

# keep only samples with non-NA cancer_group and cohort
# add cancer_group_cohort column
td_htl_dfs <- lapply(td_htl_dfs, function(x) {
  # assert all Kids_First_Biospecimen_IDs are unique
  stopifnot(identical(nrow(x), length(unique(x$Kids_First_Biospecimen_ID))))

  fx <- filter(x, !is.na(cancer_group), !is.na(cohort))
  fx <- mutate(
    fx, cancer_group_cohort = paste(cancer_group, cohort, sep = '___'))
  return(fx)
})



# Subset tumor samples and used columns in MAF table ---------------------------
tumor_kfbids <- htl_df %>%
  filter(sample_type == 'Tumor') %>%
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
suppressMessages(gc(reset = TRUE))



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
maf_df <- maf_df %>%
  left_join(mg_qres_df, by = 'Gene')

suppressMessages(gc(reset = TRUE))



# Compute mutation frequencies -------------------------------------------------
message('Compute mutation frequencies...')
cancer_group_cohort_summary_df <- get_cg_cs_tbl(td_htl_dfs$overall_htl_df)

# nf = n_samples filtered
nf_cancer_group_cohort_summary_df <- cancer_group_cohort_summary_df %>%
  filter(n_samples >= 5)

mut_freq_tbl_list <- lapply(
  1:nrow(nf_cancer_group_cohort_summary_df),
  function(i) {
    cgcs_row <- nf_cancer_group_cohort_summary_df[i, ]
    stopifnot(identical(nrow(cgcs_row), as.integer(1)))
    stopifnot(identical(length(cgcs_row$cohort_list), as.integer(1)))
    stopifnot(is.character(cgcs_row$cohort_list[[1]]))

    c_cancer_group <- cgcs_row$cancer_group
    c_cohorts <- cgcs_row$cohort_list[[1]]
    stopifnot(identical(paste(c_cohorts, collapse = '&'), cgcs_row$cohort))
    message(paste(c_cancer_group, cgcs_row$cohort))

    res_df <- get_cg_ch_mut_freq_tbl(
      maf_df, td_htl_dfs$overall_htl_df, td_htl_dfs$primary_htl_df,
      td_htl_dfs$relapse_htl_df, c_cancer_group, c_cohorts)
    return(res_df)
  }
)

m_mut_freq_tbl <- bind_rows(mut_freq_tbl_list) %>%
  distinct()

# assert all rows are unique
stopifnot(identical(nrow(m_mut_freq_tbl),
                    sum(sapply(mut_freq_tbl_list, nrow))))

write_tsv(m_mut_freq_tbl,
          file.path(tables_dir, 'snv-consensus-annotated-mut-freq.tsv'))
