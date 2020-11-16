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
# was used to subset clinical file
lgat_specimens_df<-read_tsv(file.path(subset_dir, "lgat_metadata.tsv"))

# Filter to RNA-Seq samples
lgat_wgs_df <- lgat_specimens_df %>%
  dplyr::filter(experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID)

# get putative oncogene fusion list
putativeFusion <- readr::read_tsv(file.path(root_dir,
                                            "analyses",
                                            "fusion-summary",
                                            "results",
                                            "fusion_summary_lgat_foi.tsv")) %>%
  # filter for LGAT RNA-Seq biospecimen
  dplyr::filter(Kids_First_Biospecimen_ID %in% lgat_wgs_df$Kids_First_Biospecimen_ID)

# fusion list of interest
fusionOI <- jsonlite::fromJSON(file.path(root_dir,
                                         "analyses",
                                         "molecular-subtyping-LGAT",
                                         "input",
                                         "fusionOI_list.json"))

# collapse gene list with "|" for easier grep
MAPK_fused_gene <- paste(fusionOI$MAPK$gene,collapse = "|")
RTK_fused_gene <- paste(fusionOI$RTK$gene,collapse = "|")
FGFR_fused_gene <- paste(fusionOI$FGFR$gene,collapse = "|")
MYB_fused_gene <- paste(fusionOI$MYB$gene,collapse = "|")

subsetFusion <- putativeFusion %>%
  mutate(
    # get KIAA1549--BRAF status
    KIAA_BRAF_fus = case_when(
      # canonical BRAF fusion 
      grepl("1",.$`KIAA1549--BRAF`) ~ "Yes",
      TRUE ~ "No"
    ),
    # other MAPK fusion status
    MAPK_fus = case_when(
      # Fusion with BRAF and RAF1 
      rowSums(dplyr::select(putativeFusion,dplyr::matches(MAPK_fused_gene))) > 0 &
        # remove biospecimens with canonical fusion
        # they are a separate subtype as shown above
        !grepl("1",.$`KIAA1549--BRAF`) ~ "Yes",
      TRUE ~ "No"
    ),
    # RTK fusion status
    RTK_fus = case_when(
      # fusion in any RTK gene ALK|ROS1|NTRK2|NTRK1|NTRK3|PDGFRA
      rowSums(dplyr::select(putativeFusion,dplyr::matches(RTK_fused_gene))) > 0 ~ "Yes",
      TRUE ~ "No"
    ),
    # get FGFR fusion status
    FGFR_fus = case_when(
      # fusion in any FGFR gene FGFR1|FGFR2
      rowSums(dplyr::select(putativeFusion,dplyr::matches(FGFR_fused_gene))) > 0 ~ "Yes",
      TRUE ~ "No"
    ),
    # get MYB mutation status
    MYB_fus = case_when(
      # fusion in any MYB gene MYB|MYBL1
      rowSums(dplyr::select(putativeFusion,dplyr::contains(MYB_fused_gene))) > 0 ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                KIAA_BRAF_fus,
                MAPK_fus,
                RTK_fus,
                FGFR_fus,
                MYB_fus)

# save to subset folder
write_tsv(subsetFusion,file.path(subset_dir, "LGAT_fusion_subset.tsv"))

