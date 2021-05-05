# Take in oncoprint-goi-lists-OpenPBTA.csv and create a goi file for each
# column with associated genes of interest for each specified broad histology
#
#
# Chante Bethell for CCDL 2021
#
# # #### USAGE
# This script is intended to be sourced in the script as follows:
# 
# Rscript --vanilla 00-prepare-goi-lists.R \


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
all_goi_df <- readr::read_csv(file.path("data", 
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
