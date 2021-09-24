# Purpose: Generate tables of independent rna-seq specimens 

# load libraries
library(magrittr)
library(dplyr)

# base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "independent-samples")
out_dir <- file.path(analysis_dir, "results")
dir.create(out_dir, showWarnings = F, recursive = T)

# source function
source(file.path(analysis_dir, "util", "independent_rna_samples.R"))

# read histology file
histology_df <- readr::read_tsv(file.path(root_dir, 'data/histologies.tsv'))

# randomize rows of histology file to avoid selection bias
set.seed(100)
histology_df <- histology_df[sample(nrow(histology_df)), ]

# Read in DNA independent sample list to match to rna samples
# So that independent RNA samples match the DNA samples
independent_dna_sample_df_each <- readr::read_tsv("results/independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv")
independent_dna_sample_df_all <- readr::read_tsv("results/independent-specimens.wgswxspanel.primary-plus.tsv")

# Filter to only samples from tumors, where composition is known to be Solid Tissue or Bone Marrow
# all RNA samples
histology_df <- histology_df %>%
  dplyr::filter(sample_type == "Tumor", 
                composition == "Solid Tissue" | composition == "Bone Marrow") 

# write independent sample outputs for independent level of each cohort 
rnaseq_primary_each_file <- file.path(out_dir, "independent-specimens.rnaseq.primary.eachcohort.tsv")
independent_rna_primary_each <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_each,
                          independent_level = "each-cohort",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary", 
                          seed = 2020) %>% 
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primary_each_file)

rnaseq_relapse_each_file <- file.path(out_dir, "independent-specimens.rnaseq.relapse.eachcohort.tsv")
independent_rna_relapse_each <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_each,
                          independent_level = "each-cohort",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "relapse", 
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_relapse_each_file)

rnaseq_primplus_each_file <- file.path(out_dir, "independent-specimens.rnaseq.primary-plus.eachcohort.tsv")
independent_rna_primary_plus_each <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_each,
                          independent_level = "each-cohort",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary_plus",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primplus_each_file)

# write independent sample outputs for independent level of all cohorts 
rnaseq_primary_all_file <- file.path(out_dir, "independent-specimens.rnaseq.primary.tsv")
independent_rna_primary_all <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primary_all_file)

rnaseq_relapse_all_file <- file.path(out_dir, "independent-specimens.rnaseq.relapse.tsv")
independent_rna_relapse_all <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "relapse",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_relapse_all_file)

rnaseq_primplus_all_file <- file.path(out_dir, "independent-specimens.rnaseq.primary-plus.tsv")
independent_rna_primary_plus_all <- histology_df %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary_plus",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primplus_all_file)
