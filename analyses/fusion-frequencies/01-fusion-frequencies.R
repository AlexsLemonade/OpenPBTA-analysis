# Read filtered fusion calls to convert into json format for OpenTarget portal
# attached gene and disease annotation and gatheres frequency at FusionName or
# Gene_Symbol level

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--fusion_file"), type = "character",
              help = "Input filtered fusion from `fusion_filtering` module "),
  make_option(c("-a", "--alt_id"), type = "character",
              help = "columnname from fusion_file to be used as Alt_ID or multiple comma separated columnames"),
  make_option(c("-c", "--input_histologies"), type  = "character",
              help = "Input histologies.tsv "),
  make_option(c("-p", "--primary_independence_list"), type = "character",
              help = "Input independence list for primary samples"),
  make_option(c("-r", "--relapse_independence_list"), type = "character",
              help = "Input independence list for relapse samples"),
  make_option(c("-o", "--output_filename"), type = "character", 
              default = "putative-oncogene-fusion-freq",
              help = "Output filename suffix ")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
fusion_file <- opt$fusion_file
input_histologies <- opt$input_histologies
primary_independence_list <- opt$primary_independence_list
relapse_independence_list <- opt$relapse_independence_list
output_filename <- opt$output_filename
alt_id <- unlist(strsplit(opt$alt_id,","))

# Create output dir ------------------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results, plots, and subset files directories
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "fusion-frequencies")
results_dir <- file.path(module_dir, "results")

# Source function to gather counts and frequency calculation per cohort and cancer_group
source(file.path(module_dir,"utils","freq_counts.R"))

# Create results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read fusion, independent sample list, and histology data ------------------------
message('Read data...')
htl_df <- read_tsv(input_histologies, guess_max = 100000,
                   col_types = cols(.default = col_guess()))
# assert no Kids_First_Biospecimen_ID or Kids_First_Participant_ID is NA
stopifnot(identical(
  sum(is.na(select(htl_df, Kids_First_Biospecimen_ID,
                   Kids_First_Participant_ID))),
  as.integer(0)))

fusion_df <- read_tsv(fusion_file)
# assert all records have Sample
stopifnot(identical(sum(is.na(fusion_df$Sample)), as.integer(0)))

# primary independent sample data frame
primary_indp_sdf <- read_tsv(primary_independence_list,
  col_types = cols(
    .default = col_guess()))

# relapse independent samples
relapse_indp_sdf <- read_tsv(relapse_independence_list,
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
  return(fx)
})

# Subset tumor samples to filter alteration df ---------------------------
tumor_kfbids <- htl_df %>%
  filter(sample_type == 'Tumor') %>%
  pull(Kids_First_Biospecimen_ID)


fusion_df <- fusion_df  %>%
  # We want to keep track of the gene symbols for each sample-fusion pair
  dplyr::select(Sample,
                Kids_First_Participant_ID, 
                FusionName, 
                Gene1A, 
                Gene1B, 
                Gene2A, 
                Gene2B, 
                Fusion_Type,
                DomainRetainedGene1A,
                DomainRetainedGene1B,
                reciprocal_exists,
                annots,
                BreakpointLocation,
                ends_with("anno")) %>%
  # Update reciprocal_exits to fusion with atleast 1 kinase gene involved
  mutate(
         reciprocal_exists_kinase = 
           # if either gene are kinase we will have "Yes" or "No" values
           # in DomainRetainedGene1A and DomainRetainedGene1B
           if_else(is.na(DomainRetainedGene1A) &
                     is.na(DomainRetainedGene1B) ,
                   # reciprocl_exists is a logical column
                   FALSE,
                   reciprocal_exists)) %>%
  dplyr::rename("Kids_First_Biospecimen_ID"="Sample",
         "Kinase_domain_retained_Gene1A" = "DomainRetainedGene1A",
         "Kinase_domain_retained_Gene1B" = "DomainRetainedGene1B",
         "Reciprocal_exists_either_gene_kinase" = "reciprocal_exists_kinase")%>%
  # replace NA to "" in columns that have NA
  replace_na(list("Kinase_domain_retained_Gene1A"="",
                  "Kinase_domain_retained_Gene1B"="",
                  "Reciprocal_exists_either_gene_kinase"="",
                  "Gene1A_anno"="",
                  "Gene1B_anno"="",
                  "Gene2A_anno"="",
                  "Gene2B_anno"="",
                  "Fusion_anno"=""))%>%
  # We want a single column that contains the gene symbols
  tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                key = gene_position, value = gene_symbol) %>%
  filter(!is.na(gene_symbol)) %>%
  dplyr::rename(Gene_Position = gene_position,
                Gene_Symbol = gene_symbol)

print(alt_id)

if (opt$alt_id == "FusionName,Fusion_Type") {
  fusion_df <- fusion_df %>%
    mutate( Alt_ID=paste(!!as.name(alt_id),sep="_"))
  keep_cols <- c("FusionName","Fusion_Type","Gene_Symbol", "Gene_Position",
                 "Fusion_anno","BreakpointLocation", "annots",
                 "Kinase_domain_retained_Gene1A",
                 "Kinase_domain_retained_Gene1B",
                 "Reciprocal_exists_either_gene_kinase",
                 "Gene1A_anno","Gene1B_anno","Gene2A_anno","Gene2B_anno")
  
} else if (opt$alt_id == "Gene_Symbol") {
  fusion_df <- fusion_df %>%
    mutate( Alt_ID = Gene_Symbol )
  keep_cols <-c("Gene_Symbol")

}

rm(tumor_kfbids)


# Compute mutation frequencies -------------------------------------------------
message('Compute mutation frequencies...')
cancer_group_cohort_summary_df <- get_cg_cs_tbl(td_htl_dfs$overall_htl_df)

# nf = n_samples filtered
nf_cancer_group_cohort_summary_df <- cancer_group_cohort_summary_df %>%
  filter(n_samples >= 5)

# Run get_cg_ch_mut_freq_tbl() per cancer_group and cohort
fus_freq_tbl_list <- lapply(
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
    
    # Call function for each cohort and cancer_group
    res_df <- get_cg_ch_mut_freq_tbl(
      fusion_df, td_htl_dfs$overall_htl_df, td_htl_dfs$primary_htl_df,
      td_htl_dfs$relapse_htl_df, c_cancer_group, c_cohorts) %>%
      # if Alt_ID has no matched Patients_in_dataset
      filter(!is.na(Patients_in_dataset))
    
    if ( nrow(res_df) != 0 ){
      # merge fusion dataframe with the counts and frequencies in each cancer_group
      # - per cohort 
      # - per primary samples in cohort
      # - per relapse samples in cohort
      fusion_df <- fusion_df %>% 
        select(-Kids_First_Biospecimen_ID,-Kids_First_Participant_ID) %>%
        unique() %>%
        inner_join(res_df,by="Alt_ID") 

      return(fusion_df)
    }
  }
)

m_fus_freq_tbl <- bind_rows(fus_freq_tbl_list) %>%
  distinct()

m_fus_freq_tbl <- m_fus_freq_tbl %>%
  mutate(Dataset = if_else(Dataset == 'PBTA&GMKF&TARGET',
                           true = 'PedOT', false = Dataset))

stopifnot(identical(sum(is.na(m_fus_freq_tbl)), as.integer(0)))

m_fus_freq_tbl <- m_fus_freq_tbl %>%
  select(keep_cols, Disease,
         Total_alterations_Over_Patients_in_dataset,
         Frequency_in_overall_dataset,
         Total_primary_tumors_mutated_Over_Primary_tumors_in_dataset,
         Frequency_in_primary_tumors,
         Total_relapse_tumors_mutated_Over_Relapse_tumors_in_dataset,
         Frequency_in_relapse_tumors) %>%
  unique() %>%
  write_tsv(file.path(results_dir, paste0(output_filename,'.tsv')))


