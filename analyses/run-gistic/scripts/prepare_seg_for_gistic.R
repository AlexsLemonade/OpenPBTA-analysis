# prepare_seg_for_gistic.R
#
# 
# Purpose: Generate seg files that are compatible with GISTIC 
# Currently, if we use the cnv consensus seg file as is, since we define copy.num=NA for neutral calls
# a lot of NAs are present in the file and will cause GISTIC to error out. 
# Those changes were made for the focal CN module so in order not to change that but also fixed the 
# issue in GISTIC, we generated a consensus seg file just for GISTIC by using `tumor_ploidy` as `copy.num`


# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Parse command line options
library(optparse)
library(tidyverse)

option_list <- list(
  make_option(
    c("-i", "--in_consensus"),
    type = "character",
    default = NULL,
    help = "cnv consensus seg file used for focal CN (.seg)",
  ),
  make_option(
    c("-o", "--out_consensus"),
    type = "character",
    default = NULL,
    help = "cnv consensus seg file used for GISTIC (.seg)"
  ),
  make_option(
    c("-c", "--histology"),
    type = "character",
    default = NULL, 
    help = "histology file to match the ploidy of samples (.TSV)"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))


# Read the cnv consensus file with copy.num=NA for neutral calls
cnv_consensus <- readr::read_tsv(opts$in_consensus)
# Read in histology file for ploiidy information
metadata <- readr::read_tsv(opts$histology, guess_max=100000)

# Use the tumor ploidy for 
cnv_consensus_for_gistic <- metadata %>% 
  dplyr::select(Kids_First_Biospecimen_ID, tumor_ploidy) %>%
  dplyr::rename(ID = Kids_First_Biospecimen_ID) %>%
  dplyr::inner_join(cnv_consensus, by="ID") %>%
  dplyr::mutate(copy.num = dplyr::case_when(is.na(copy.num) ~ tumor_ploidy, TRUE ~ copy.num)) %>%
  dplyr::select(-tumor_ploidy)

# write out the file for GISTIC usage
readr::write_tsv(cnv_consensus_for_gistic,opts$out_consensus)

