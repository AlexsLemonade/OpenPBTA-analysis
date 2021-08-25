suppressPackageStartupMessages({
  library(tidyverse)
})

# Comparing EXTEND Scores with histology molecular subtypes 
# base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 
analysis_dir <- file.path(root_dir, "analyses", "telomerase-activity-prediction")

# source function
source(file.path(analysis_dir, "util", "plot_by_molecular_subtype.R"))

# histology file
hist_file <- file.path(root_dir, "data", "pbta-histologies.tsv") 

# telomerase activity-prediction for Stranded_FPKM data
telomerase_scores <- file.path(analysis_dir, "results", "TelomeraseScores_PTBAStranded_FPKM.txt") 

# output file
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

## Reading the clinical data
hist_file <- read_tsv(hist_file, guess_max = 10000)    
hist_file <- hist_file %>%
  rename("SampleID" = "Kids_First_Biospecimen_ID")

## Reading Stranded FPKM telomerase scores
telomerase_scores <- read.delim(telomerase_scores)  

### Merging Clinical data with the Telomerase scores
telomerase_scores <- hist_file %>%
  select(SampleID, short_histology, broad_histology, molecular_subtype) %>%
  inner_join(telomerase_scores, by = "SampleID")

# filter NA, to be classified and reduce to sample count >= 3
telomerase_scores <- telomerase_scores %>%
  filter(!is.na(molecular_subtype),
         !grepl("To be classified", molecular_subtype)) %>%
  group_by(broad_histology, molecular_subtype) %>%
  mutate(n_samples = n()) %>%
  filter(n_samples >= 3)

# filter NA, to be classified and molecular subtype count < 2 
subtypes <- telomerase_scores %>%
  select(broad_histology, molecular_subtype) %>%
  unique() %>%
  group_by(broad_histology) %>%
  mutate(n_subtypes = n()) %>%
  filter(n_subtypes > 1)
telomerase_scores <- telomerase_scores %>%
  filter(molecular_subtype %in% subtypes$molecular_subtype)

# apply function
# create plot of NormEXTENDScores per molecular subtype per histology 
plyr::d_ply(telomerase_scores, 
            .variables = "broad_histology", 
            .fun = function(x) plot_by_molecular_subtype(scores_mat = x, plots_dir, results_dir))

