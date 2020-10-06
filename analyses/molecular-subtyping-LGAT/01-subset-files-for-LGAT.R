# K S Gaonkar 2020
# Subsets the consensus mutation file to biospecimens where short_histology == LGAT

library(tidyverse)

# Look for git root folder
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# get subset folder
subset_dir <- file.path(root_dir, 
                        "analyses",
                        "molecular-subtyping-LGAT",
                        "lgat-subset")

# create if doesn't exist
if (!dir.exists(subset_dir)) {
  dir.create(subset_dir)
}

# File from 00-LGAT-select-pathology-dx that is used for the pathology diagnosis
# inclusion/exclusion criteria
path_dx_list <- jsonlite::fromJSON(
  file.path(subset_dir, 
            "lgat_subtyping_path_dx_strings.json")
)

# clinical file
clinical <- read_tsv(file.path(root_dir, 
                              "data",
                              "pbta-histologies.tsv"), 
                    guess_max = 10000)
# consensus mutation data
consensusMutation <- read_tsv(file.path(root_dir,
                                        "data",
                                        "pbta-snv-consensus-mutation.maf.tsv.gz"))

# Filter to tumor samples that should be included on the basis of pathology
# diagnosis
lgat_specimens_df <- clinical %>%
  filter(str_detect(str_to_lower(pathology_diagnosis),  # Inclusion criteria
                    paste0(path_dx_list$include_path_dx, collapse = "|")),
         # Exclusion criteria
         str_detect(str_to_lower(pathology_diagnosis),
                    paste0(path_dx_list$exclude_path_dx, collapse = "|"),
                    negate = TRUE),
         # Tumors
         sample_type == "Tumor",
         composition == "Solid Tissue")

# Write this intermediate file to the subset directory as it allows for
# inspection
write_tsv(lgat_specimens_df, file.path(subset_dir, "lgat_metadata.tsv"))

# Filter to WGS samples
lgat_wgs_df <- lgat_specimens_df %>%
  filter(experimental_strategy == "WGS") %>%
  select(Kids_First_Biospecimen_ID)

# Filter consensus mutation files for LGAT subset
consensusMutationSubset <- consensusMutation %>%
  # find lgat samples
  filter(Tumor_Sample_Barcode %in% lgat_wgs_df$Kids_First_Biospecimen_ID) %>%
  # look for BRAF V600E mutations and make a column for BRAF_V600E
  filter(Hugo_Symbol == "BRAF" & HGVSp_Short == "p.V600E") %>%
  # select tumor sample barcode
  select(Tumor_Sample_Barcode, HGVSp_Short) %>% 
  distinct() %>%
  # join other WGS LGAT samples
  full_join(lgat_wgs_df,
            by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  # get BRAF_V600E status
  mutate(BRAF_V600E= case_when(
    HGVSp_Short == "p.V600E" ~ "Yes",
    is.na(HGVSp_Short) ~ "No"
  ))

# remove consensusMutation
rm(consensusMutation)

# save to subset folder
write_tsv(consensusMutationSubset,file.path(subset_dir, "LGAT_snv_subset.tsv"))
