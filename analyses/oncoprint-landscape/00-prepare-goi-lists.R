# Take in oncoprint-goi-lists-OpenPBTA.csv and create a goi file for each
# column with associated genes of interest for each specified broad histology.
# Also creates a table for mapping between cancer_group and the appropriate GOI
# list.
#
#
# Chante Bethell for CCDL 2021
#
# USAGE:
# 
# Rscript --vanilla 00-prepare-goi-lists.R


#### Set Up --------------------------------------------------------------------

library(dplyr)
library(stringr)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to directory that contains the manually curated genes of interest list
data_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "data")

# Each histology has it's own column, and there are source columns
all_goi_df <- readr::read_csv(file.path(data_dir, 
                                        "oncoprint-goi-lists-OpenPBTA.csv"))

# Drop the source columns
all_goi_df <- all_goi_df %>% 
  select(-contains("Source"))

# Now each column will be a broad histology
for (col_iter in 1:ncol(all_goi_df)) {
  
  # The broad histology is the column name, but let's make it all lowercase
  # and replace spaces with hyphens for use as part of the output file name
  broad_histology <- str_to_lower(str_replace_all(
      colnames(all_goi_df)[col_iter],
      pattern = " ",
      replacement = "-"
    ))
  
  # Create the output file name
  output_file <- file.path(data_dir, str_c(broad_histology, "_goi_list.tsv"))
  
  # Write the current column to file, removing any NA values
  all_goi_df[, col_iter] %>%
    tidyr::drop_na() %>%
    distinct() %>%
    readr::write_tsv(output_file)

}

#### For specifically getting cancer group-specific counts of the alterations --
# We will something that maps between cancer group and gene of interest list.
# We'll also include the oncoprint_group from the palette to make inspection
# easier.

palette_df <- readr::read_tsv(file.path(root_dir, 
                                        "figures", 
                                        "palettes", 
                                        "broad_histology_cancer_group_palette.tsv"))

cg_df <- palette_df %>% 
  select(oncoprint_group, cancer_group, cancer_group_display) %>% 
  filter(!is.na(oncoprint_group)) %>%
  mutate(goi_list_file = case_when(
    oncoprint_group == "Other CNS" ~ "data/other_goi_list.tsv",
    oncoprint_group == "Low-grade astrocytic tumor" ~ "data/lgat_goi_list.tsv",
    oncoprint_group == "Embryonal tumor" ~ "data/embryonal-tumor_goi_list.tsv",
    oncoprint_group == "Diffuse astrocytic and oligodendroglial tumor" ~ "data/hgat_goi_list.tsv",
    TRUE ~ NA_character_
  ))

readr::write_tsv(cg_df, 
                 file.path(data_dir, "cancer_group_goi_list_mapping.tsv"))
