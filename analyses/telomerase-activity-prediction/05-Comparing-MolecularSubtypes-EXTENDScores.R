suppressPackageStartupMessages({
  library(tidyverse)
})

# Comparing EXTEND Scores with histology molecular subtypes 
# base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 
analysis_dir <- file.path(root_dir, "analyses", "telomerase-activity-prediction")

# source function
source(file.path(analysis_dir, "util", "boxplot_by_molecular_subtype.R"))

# histology file
hist_file <- file.path(root_dir, "data", "pbta-histologies.tsv") 

# telomerase activity-prediction for Stranded_FPKM data
telomerase_scores <- file.path(analysis_dir, "results", "TelomeraseScores_PTBAStranded_FPKM.txt") 

# output file
output_file <- file.path(analysis_dir, "plots", "PBTA_MolecularSubtypes.pdf")

## Reading the clinical data
hist_file <- readr::read_tsv(hist_file, guess_max = 10000)    
hist_file <- hist_file %>%
  dplyr::rename("SampleID" = "Kids_First_Biospecimen_ID")

## Reading Stranded FPKM telomerase scores
telomerase_scores <- read.delim(telomerase_scores)  

### Merging Clinical data with the Telomerase scores
telomerase_scores <- hist_file %>%
  dplyr::select(SampleID, short_histology, broad_histology, molecular_subtype) %>%
  inner_join(telomerase_scores, by = "SampleID")

# filter NA, to be classified and molecular subtype count < 2 (required for p-value calculations)
telomerase_scores <- telomerase_scores %>%
  filter(!is.na(molecular_subtype),
         !grepl("To be classified", molecular_subtype)) %>%
  group_by(molecular_subtype) %>%
  mutate(n = n()) %>%
  filter(n > 1)

# apply function
# create boxplot of NormEXTENDScores per molecular subtype per histology 
pdf(output_file, width = 7, height = 5)
plyr::d_ply(telomerase_scores, 
            .variables = "short_histology", 
            .fun = function(x) boxplot_by_molecular_subtype(scores_mat = x))
dev.off()

