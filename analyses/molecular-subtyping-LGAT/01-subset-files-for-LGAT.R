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

# combined list for SNV of interest
snvOI <- jsonlite::fromJSON(file.path(root_dir,
                                      "analyses",
                                      "molecular-subtyping-LGAT",
                                      "input",
                                      "snvOI_list.json"))

# collase multiple hotspots in genes with "|" so easy grep calls
BRAF_hotspot <-paste(snvOI$BRAF_V600E$hotspot[!is.na( snvOI$BRAF_V600E$hotspot)],collapse = "|")
FGFR_hotspot <-paste(snvOI$FGFR$hotspot[!is.na( snvOI$FGFR$hotspot)],collapse = "|")
IDH_hotspot<- paste(snvOI$IDH$hotspot[!is.na( snvOI$IDH$hotspot)],collapse = "|")


# Filter to tumor samples that should be included on the basis of pathology
# diagnosis
lgat_specimens_df <- clinical %>%
  dplyr::filter(str_detect(str_to_lower(pathology_diagnosis),  # Inclusion criteria
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
  dplyr::filter(experimental_strategy == "WGS") %>%
  dplyr::select(Kids_First_Biospecimen_ID)

# Filter consensus mutation files for LGAT subset
consensusMutationSubset <- consensusMutation %>%
  # find lgat samples
  dplyr::filter(Tumor_Sample_Barcode %in% lgat_wgs_df$Kids_First_Biospecimen_ID) %>%
  # select tumor sample barcode
  dplyr::select(Tumor_Sample_Barcode,
                Hugo_Symbol,
                HGVSp_Short,
                DOMAINS,
                Variant_Classification) %>% 
  dplyr::distinct() %>%
  # join other WGS LGAT samples
  dplyr::full_join(lgat_wgs_df,
            by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  mutate(
    # get BRAF mutation status
    BRAF_V600E_mut = case_when(
      # canonical mutations V600E
      HGVSp_Short %in% snvOI$BRAF_V600E$canonical[!is.na(snvOI$BRAF_V600E$canonical)] &
        Hugo_Symbol=="BRAF" ~ "Yes",
      # hotspot mutations in p.600 and p.599
      grepl(BRAF_hotspot,HGVSp_Short) &
        Hugo_Symbol=="BRAF" &
        Variant_Classification != "Silent" ~ "Yes",
      # and kinase domain mutation for non-canonical mutation 
      grepl("PF07714",DOMAINS) & 
        Hugo_Symbol=="BRAF" & 
        Variant_Classification != "Silent" ~ "Yes",
      TRUE ~ "No"
    ),
    # get NF1 mutation status
    NF1_mut = case_when(
      # all mutations in MAPK genes
      Hugo_Symbol %in% snvOI$NF1$gene & 
        Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation") ~ "Yes",
      TRUE ~ "No"
    ),
    # other MAPK mutation status
    MAPK_mut = case_when(
      # all mutations in MAPK genes
      Hugo_Symbol %in% snvOI$MAPK$gene ~ "Yes",
      TRUE ~ "No"
    ),
    # RTK mutation status
    RTK_mut = case_when(
      # all mutations in RTK genes
      Hugo_Symbol %in% snvOI$RTK$gene ~ "Yes",
      TRUE ~ "No"
    ),
    # get FGFR mutation status
    FGFR_mut = case_when(
      # canonical mutations
      HGVSp_Short %in% snvOI$FGFR$canonical[!is.na(snvOI$FGFR$canonical)] &
        Hugo_Symbol=="FGFR1"~ "Yes",
      # hotspot mutations 
      grepl(FGFR_hotspot,HGVSp_Short) &
        Hugo_Symbol=="FGFR1" ~ "Yes",
      TRUE ~ "No"
    ),
    # get IDH mutation status
    IDH_mut = case_when(
      # hostspot mutations
      grepl(IDH_hotspot,HGVSp_Short) & 
        Hugo_Symbol %in% snvOI$IDH$gene ~ "Yes",
      TRUE ~ "No"
    ),
    # get H3F3A mutation status
    H3F3A_mut = case_when(
      # canonical mutations
      HGVSp_Short %in% snvOI$H3F3A$canonical & Hugo_Symbol %in% "H3F3A" ~ "Yes",
      TRUE ~ "No"
    ),
    # get HIST1H3B mutation status
    HIST1H3B_mut = case_when(
      # canonical mutations
      HGVSp_Short %in% snvOI$HIST1H3B$canonical & Hugo_Symbol %in% "HIST1H3B" ~ "Yes",
      TRUE ~ "No"
    ),
    # get HIST1H3C mutation status
    HIST1H3C_mut = case_when(
      # canonical mutations
      HGVSp_Short %in% snvOI$HIST1H3C$canonical & Hugo_Symbol %in% "HIST1H3C" ~ "Yes",
      TRUE ~ "No"
    ),
  ) %>%
  dplyr::select(Tumor_Sample_Barcode,
                BRAF_V600E_mut,
                FGFR_mut,
                IDH_mut,
                H3F3A_mut,
                HIST1H3B_mut,
                HIST1H3C_mut)

# remove consensusMutation
rm(consensusMutation)

# save to subset folder
write_tsv(consensusMutationSubset,file.path(subset_dir, "LGAT_snv_subset.tsv"))

