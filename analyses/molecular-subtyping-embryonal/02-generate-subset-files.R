# Stephanie J. Spielman and Jaclyn Taroni for ALSF CCDL 2020
#
# This script subsets the files required for subtyping non-MB and non-ATRT
# embryonal tumors. A sample will be included if either of the following
# conditions are met: 1) the sample is labeled as an embryonal tumor
# (broad_histology) but NOT an MB or ATRT tumor (disease_type_old) OR 2)
# the sample contains a TTYH1 fusion (5' partner)

library(tidyverse)

#### Directories ---------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-embryonal")

# Set path to directory where we will save the subset files
subset_dir <- file.path(module_dir, "subset-files")
if (!dir.exists(subset_dir)) {
  dir.create(subset_dir)
}

# Results directory
results_dir <- file.path(module_dir, "results")

#### Read in files -------------------------------------------------------------

# File that contains the relevant biospecimen identifiers -- this was generated
# in 01-samples-to-subset
subset_id <-
  read_tsv(file.path(results_dir,
                     "biospecimen_ids_embryonal_subtyping.tsv")) %>%
  pull(Kids_First_Biospecimen_ID)


