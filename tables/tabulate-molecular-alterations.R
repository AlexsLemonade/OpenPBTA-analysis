# S. Spielman and J. Taroni 2023 for CCDL
# 
# This script tabulates, for _all samples_ in OpenPBTA cohort, the specific alterations present in 
# each sample for genes that are present in the manuscript's oncoprint plots.
# Aspects of this script are inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R


# Define pipe
`%>%` <- dplyr::`%>%`

# We specifically are interested in showing the Figure 2 (and S3B which has same genes as 2D) alterations but for all tumors.
# The gene_list below contains the top GOIs identified in `scripts/figures/fig2_figS3-oncoprint-landscape.R`, sorted and de-duplicated
# For further reference, see Figure 2:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/aa929753fca0019294571ec813e6ae7224b1d8b8/figures/pngs/figure2.png
gene_list <- c("ACVR1", "APC", "ATM", "ATRX", "BCOR", "BRAF", 
               "C11orf95", "CDK6", "CTDNEP1", "CTNNB1", "DDX3X", 
               "EGFR", "ERBB4", "EWSR1", "FGFR1", "FGFR2", "FOXO3", 
               "H3F3A", "IGF1R", "KDM6A", "KIAA1549", "KIT", 
               "KMT2C", "KMT2D", "KRAS", "KSR2", "MAML2", "MAMLD1", 
               "MET", "MIR512-2", "MYB", "MYCN", "NF1", "NF2", 
               "NTRK2", "NTRK3", "PDGFRA", "PIK3CA", "PIK3R1", 
               "PPM1D", "PRKAR1A", "PTCH1", "PTEN", "PTPN11", 
               "QKI", "RAF1", "RELA", "ROS1", "SETD2", "SMARCA4", 
               "TACC1", "TCF4", "TERT", "TP53", "YAP1", "ZIC1")

# Later we'll need to use a list to fill in NAs with tidyr::replace_na()
gene_fill_list <- list()
for (gene in gene_list) gene_fill_list[[gene]] <-  "None"

## Set up paths ----------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
palette_file <- file.path(root_dir, 
                          "figures", 
                          "palettes", 
                          "broad_histology_cancer_group_palette.tsv")
metadata_file <- file.path(data_dir, "pbta-histologies.tsv")
maf_file <- file.path(data_dir, 
                      "pbta-snv-consensus-mutation.maf.tsv.gz")
hotspots_maf_file <- file.path(data_dir, 
                               "pbta-snv-scavenged-hotspots.maf.tsv.gz")
cnv_autosomes_file <- file.path(data_dir, 
                                "consensus_seg_annotated_cn_autosomes.tsv.gz")
cnv_xy_file <- file.path(data_dir, 
                         "consensus_seg_annotated_cn_x_and_y.tsv.gz")
fusion_file <- file.path(data_dir, 
                         "pbta-fusion-putative-oncogenic.tsv")

# Output CSV file
zenodo_csv_file <- file.path(root_dir, 
                             "tables", 
                             "zenodo-upload",
                             "openpbta-molecular-alterations.tsv")

## Read in data ----------------------------------------------------------------

# metadata
histologies_df <- readr::read_tsv(metadata_file, guess_max = 10000)
palette_df <- readr::read_tsv(palette_file) 

# maf and hotspots data
maf_df <- readr::read_tsv(maf_file)
hotspot_df <- readr::read_tsv(hotspots_maf_file) 

# fusion data
fusion_df <- readr::read_tsv(fusion_file)

# CNV data: autosomes and XY
cnv_df <- readr::read_tsv(cnv_autosomes_file) %>%
  dplyr::bind_rows(
    readr::read_tsv(cnv_xy_file)
  )

## Prepare the maf data --------------------------------------------------------

# Only keep maf columns relevant for annotating specific alterations
maf_keep_cols <-c("Tumor_Sample_Barcode",
                  "Hugo_Symbol", 
                  "Variant_Type",
                  "HGVSp", 
                  "HGVSc")

# subset columns first to enable combining dfs; they have some different 
#  data types in other columns we're getting rid of
hotspot_df <- hotspot_df %>%
  dplyr::select(maf_keep_cols)

# Final MAF alterations
maf_df <- maf_df %>%
  dplyr::select(maf_keep_cols) %>%
  # combine consensus MAF and hotspot MAF
  dplyr::bind_rows(hotspot_df) %>%
  # filter to only relevant genes and ids
  dplyr::filter(Hugo_Symbol %in% gene_list) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    # assign the final alteration based on which piece 
    # of information is known, in this order of priority 
    final_alteration = dplyr::case_when(
      HGVSp != "." ~ HGVSp, 
      HGVSc != "." ~ HGVSc, 
      TRUE ~ Variant_Type
    )) %>%
  # remove original alteration columns
  dplyr::select(-HGVSp, -HGVSc, -Variant_Type) %>%
  # rename `Tumor_Sample_Barcode` -> `biospecimen_id` to enable later data joins
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode)

# Add sample IDs
maf_df <- histologies_df %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                composition) %>% 
  dplyr::right_join(maf_df, 
                    by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID)


## Prepare the fusion data -----------------------------------------------------

# We'll handle fusions where reciprocal fusions exist (e.g., 
# reciprocal_exists == TRUE) separately from other fusions
fusion_reciprocal_df <- fusion_df %>%
  # Reciprocal fusions only
  dplyr::filter(reciprocal_exists) %>%
  # BSID + Gene1--Gene2
  dplyr::select(Sample, FusionName) %>%
  # filter to distinct identifier-fusion pairs
  dplyr::distinct() %>%
  dplyr::group_by(Sample, FusionName) %>%
  # Put genes in the fusions in alphabetical order to collapse
  dplyr::mutate(SortedFusionName = stringr::str_c(
    sort( stringr::str_split(FusionName, "--", simplify = TRUE) ), 
    collapse = "--")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample,
                FusionName = SortedFusionName) %>%
  dplyr::distinct()

# Combine non-reciprocal with reciprocal, and reshape into:
# Fusion_ID, Sample, Partner, Hugo_Symol
fusion_df <- fusion_df %>%
  # get non-reciprocal  rows
  dplyr::filter(!reciprocal_exists) %>%
  dplyr::select(Sample, FusionName) %>%
  # remove duplicates
  dplyr::distinct() %>%
  # combine with reciprocal rows
  dplyr::bind_rows(fusion_reciprocal_df) %>%
  # Separate the 5' and 3' genes
  tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--", 
                  # don't remove the alteration itself
                  remove = FALSE) %>%
  # Use row numbers to mark unique fusions 
  tibble::rowid_to_column("Fusion_ID") %>%
  dplyr::select(Fusion_ID, FusionName, Sample, Gene1, Gene2) %>%
  tidyr::gather("Partner", "Hugo_Symbol", c("Gene1", "Gene2")) %>%
  # remove the `Partner` and `Fusion_ID` columns that we won't need anymore:
  dplyr::select(-Partner, -Fusion_ID) %>%
  # filter to only relevant genes and ids
  dplyr::filter(Hugo_Symbol %in% gene_list) %>%
  # rename `Sample` -> `biospecimen_id` and `FusionName` -> `final_alteration`
  #  to enable later data joins
  dplyr::rename(Kids_First_Biospecimen_ID = Sample, 
                final_alteration = FusionName) 

# now can remove the reciprocal: 
rm(fusion_reciprocal_df)

# Add sample IDs
fusion_df <- histologies_df %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                composition) %>% 
  dplyr::right_join(fusion_df, 
                    by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID)

## Prepare the CNV data --------------------------------------------------------

cnv_df <- cnv_df %>%
  # filter to only relevant genes and ids
  dplyr::filter(gene_symbol %in% gene_list) %>%
  # create a single variable indicating the type of CNV we have
  dplyr::mutate(
    final_alteration = glue::glue("{cytoband}-{status}-{copy_number}")
  ) %>%
  # select only relevant columns moving forward
  dplyr::select(
    Kids_First_Biospecimen_ID = biospecimen_id, 
    Hugo_Symbol = gene_symbol, 
    final_alteration
  )

# Add sample IDs and assay type
cnv_df <- histologies_df %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                composition) %>% 
  dplyr::right_join(cnv_df, 
                    by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID)


## DNA alterations -------------------------------------------------------------

# Collapse by gene -- if more than one alteration is in a gene as detected via
# DNA assay, separate them with ;
# These will be alterations detected in *any* (at least 1) DNA assay associated 
# with this sample_id
dna_alterations_df <- maf_df %>%
  dplyr::bind_rows(cnv_df) %>%
  dplyr::group_by(sample_id, 
                  Hugo_Symbol,
                  composition) %>%
  dplyr::summarize(final_alteration = paste(sort(unique(final_alteration)),  
                                            collapse = "; ")) %>%
  dplyr::ungroup()

## RNA alterations -------------------------------------------------------------

# Collapse by gene -- if more than one fusion is the same gene, 
# separate them with ;
# These will be alterations detected in *any* RNA assay associated with this 
# sample_id
rna_alterations_df <- fusion_df %>%
  dplyr::group_by(sample_id, 
                  Hugo_Symbol,
                  composition) %>%
  dplyr::summarize(final_alteration = paste(sort(unique(final_alteration)),  
                                            collapse = "; ")) %>%
  dplyr::ungroup()

## Alterations from DNA and RNA ------------------------------------------------

alterations_df <- dna_alterations_df %>%
  # Full join of DNA and RNA alterations by sample ID and gene symbol
  dplyr::full_join(rna_alterations_df,
                   by = c("sample_id",
                          "Hugo_Symbol",
                          "composition"),
                   suffix = c("_DNA", "_RNA")) %>%
  # Combine DNA and RNA alterations into one column where applicable
  dplyr::mutate(final_alteration = dplyr::case_when(
    !is.na(final_alteration_DNA) & 
      !is.na(final_alteration_RNA) ~ paste(final_alteration_DNA,
                                           final_alteration_RNA,
                                           sep = "; "),
    is.na(final_alteration_DNA) ~ final_alteration_RNA,
    is.na(final_alteration_RNA) ~ final_alteration_DNA
  )) %>%
  # Remove the separate DNA and RNA alterations columns
  dplyr::select(sample_id,
                Hugo_Symbol,
                composition,
                final_alteration) %>%
  dplyr::distinct()

# Now spread this to get a wide version -- gene symbols become columns, and 
# genes that are not altered in that sample are filled with "None"
alterations_wide_df <- alterations_df %>%
  dplyr::group_by(sample_id, 
                  composition) %>%
  tidyr::spread(Hugo_Symbol,
                final_alteration,
                fill = "None")

## Clean up --------------------------------------------------------------------

# Remove large-ish data frames we no longer need
rm(hotspot_df, 
   cnv_df,
   maf_df,
   fusion_df,
   dna_alterations_df,
   rna_alterations_df)

## Deal with identifiers -------------------------------------------------------

# Create a data frame that tracks sample_ids, all DNA biospecimen IDs, all 
# RNA biospecimen IDs, and whether or not there is "ambiguous" mapping between
# DNA and RNA assays for a given sample_id
identifiers_df <- histologies_df %>%
  # Remove normal samples
  dplyr::filter(sample_type != "Normal") %>%
  # Get only the relevant identifiers and the experimental strategy (to be 
  # recoded)
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                composition,
                experimental_strategy) %>%
  # Code experimental strategy to more general assay type
  dplyr::mutate(assay_type = dplyr::case_when(
    experimental_strategy == "RNA-Seq" ~ "RNA",
    experimental_strategy == "WGS" ~ "DNA",
    experimental_strategy == "WXS" ~ "DNA",
    experimental_strategy == "Targeted Sequencing" ~ "DNA"
  ))

# What samples have more than one DNA or RNA assay?
multiassay_sample_ids <- identifiers_df %>%
  dplyr::group_by(sample_id, 
                  experimental_strategy,
                  composition) %>%
  dplyr::tally() %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(sample_id)

# Create a data frame where each sample_id is a row and samples with multiple
# DNA or RNA assays have the biospecimen IDs separated by ;
identifiers_df <- identifiers_df %>%
  dplyr::select(-experimental_strategy) %>%
  # Create Kids_First_Biospecimen_ID_DNA and Kids_First_Biospecimen_ID_RNA 
  # columns where multiple specimens are separated by ;
  dplyr::group_by(sample_id,
                  composition,
                  assay_type) %>%
  dplyr::summarize(biospecimen_ids = paste(sort(unique(Kids_First_Biospecimen_ID)),  
                                           collapse = "; ")) %>%
  tidyr::spread(assay_type, biospecimen_ids) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = DNA, 
                Kids_First_Biospecimen_ID_RNA = RNA) %>%
  # Add indicator of when there are multiple DNA or RNA assays per sample
  dplyr::mutate(multiple_assays_within_type = dplyr::if_else(
    sample_id %in% multiassay_sample_ids,
    true = TRUE,
    false = FALSE
  ))


## Histologies metadata to include ---------------------------------------------

metadata_df <- histologies_df %>%
  # Remove normal samples
  dplyr::filter(sample_type != "Normal") %>%
  # Include sample type, composition, germline sex estimate, broad_histology,
  # cancer group
  dplyr::select(
    sample_id,
    sample_type,
    composition,
    germline_sex_estimate,
    broad_histology,
    cancer_group
  ) %>%
  dplyr::distinct() %>%
  # Add the information about biospecimen ID
  dplyr::left_join(identifiers_df,
                    by = c("sample_id",
                           "composition")) %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID_DNA,
                Kids_First_Biospecimen_ID_RNA,
                multiple_assays_within_type,
                sample_type,
                composition,
                germline_sex_estimate,
                broad_histology,
                cancer_group,
                dplyr::everything())

# Edge cases where multiple RNA assays mean that duplicate row without germline
# sex estimate exist
retain_rows <- dplyr::if_else(
  !is.na(metadata_df$Kids_First_Biospecimen_ID_DNA) & 
    is.na(metadata_df$germline_sex_estimate), 
  FALSE, 
  TRUE)

# Remove the edge case rows
metadata_df <- metadata_df %>%
  dplyr::filter(retain_rows)


## Create and save final table -------------------------------------------------

# Add the metadata to the wide version of the alterations table
final_alterations_df <- metadata_df %>%
  dplyr::left_join(alterations_wide_df,
                   by = c("sample_id",
                          "composition")) %>%
  # If something is NA, set it to "None" with the gene fill list
  tidyr::replace_na(replace = gene_fill_list) %>%
  # Order by sample_id
  dplyr::arrange(sample_id)

# Write to file
readr::write_csv(final_alterations_df, zenodo_csv_file)

