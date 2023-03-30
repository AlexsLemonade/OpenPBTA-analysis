# S. Spielman 2023 for CCDL
# 
# This script tabulates, for _all samples_ in OpenPBTA cohort, the presence/absence of 
#  alterations depicted in the manuscript's oncoprint plots.
# This script is heavily inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R
# 

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



# Set up paths -----------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
oncoprint_dir <- file.path(root_dir,
                           "analyses", 
                           "oncoprint-landscape")

pal_file           <- file.path(root_dir, "figures", "palettes", "broad_histology_cancer_group_palette.tsv")
metadata_file      <- file.path(data_dir, "pbta-histologies.tsv")
maf_file           <- file.path(data_dir, "pbta-snv-consensus-mutation.maf.tsv.gz")
hotspots_maf_file  <- file.path(data_dir, "pbta-snv-scavenged-hotspots.maf.tsv.gz")
cnv_autosomes_file <- file.path(data_dir, "consensus_seg_annotated_cn_autosomes.tsv.gz")
cnv_xy_file        <- file.path(data_dir, "consensus_seg_annotated_cn_x_and_y.tsv.gz")
fusion_file        <- file.path(data_dir, "pbta-fusion-putative-oncogenic.tsv")

# Output CSV file
zenodo_csv_file <- file.path(root_dir, 
                             "tables", 
                             "zenodo-upload",
                             "openpbta-molecular-alterations.tsv")

#### Read in data ------------------------------------

histologies_df <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  # this incidentally also filters out "Normal" tissue from the histologies file!
  dplyr::inner_join(
    readr::read_tsv(pal_file) %>%
      dplyr::select(cancer_group, cancer_group_display, broad_histology, broad_histology_display)
  )

# Helper variable to check later steps. There are this many unique tumors + cell lines 
openpbta_total_samples <- length(unique(histologies_df$sample_id))


# Read maf and hotspots data
maf_df <- readr::read_tsv(maf_file) 
hotspots_maf_df <- readr::read_tsv(hotspots_maf_file) 

# Read fusion data
fusion_df <- readr::read_tsv(fusion_file)

# Read CNV data
cnv_df <- readr::read_tsv(cnv_autosomes_file) %>%
  dplyr::bind_rows(
    readr::read_tsv(cnv_xy_file)
  )

## Filter to relevant tumors and metadata -----------------------------------------


# TODO: UPDATE COMMENT IN THE END
# To begin, we'll update `histologies_df` to store only relevant metadata.
# As part of this, we'll note if a sample_id is ambiguous. Ambiguous sample_id's will 
# have more than 2 rows associated with it in the histologies file when looking at 
# tumor samples -- that means we won't necessarily be able to determine **when an 
# WGS/WXS assay maps to an RNA-seq assay** for the purpose of tabulating molecular alterations.
# Alternatively, ambiguous sample_id's will have multiple biospecimens of the same experimental strategy.

histologies_df <- histologies_df %>%
  # keep these columns of interest moving forward:
  dplyr::select(sample_id, 
                Kids_First_Biospecimen_ID, 
                cancer_group_display, 
                broad_histology_display, 
                experimental_strategy,
                germline_sex_estimate) %>%
  dplyr::distinct() 
  # ambiguity tracking may no longer be needed


# Check that this is equal to `openpbta_total_samples`
if (!(length(unique(histologies_df$sample_id))) == openpbta_total_samples) {
  stop("Bad sample_id parsing.")
} 

# Pull out biospecimen id column for convenience 
bs_ids <- histologies_df %>%
  dplyr::pull(Kids_First_Biospecimen_ID) %>%
  unique() 


## Prepare the maf data -------------------------------------

# Only keep maf columns relevant for annotating specific alterations
maf_keep_cols <-c("Tumor_Sample_Barcode",
                  "Hugo_Symbol", 
                  "Variant_Type",
                  "HGVSp", 
                  "HGVSc")

# Subset columns before combining dfs
maf_df <- maf_df %>%
  dplyr::select(maf_keep_cols) 
  
maf_df <- hotspots_maf_df %>%
  dplyr::select(maf_keep_cols) %>%
  # combine maf and hotspots
  dplyr::bind_rows(maf_df) %>%
  # filter to only relevant genes and ids
  dplyr::filter(Hugo_Symbol %in% gene_list,
                Tumor_Sample_Barcode %in% bs_ids) %>%
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
  dplyr::rename(biospecimen_id = Tumor_Sample_Barcode)




## Prepare the fusion data ----------------------------------

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
  dplyr::filter(Hugo_Symbol %in% gene_list,
                Sample %in% bs_ids) %>%
  # rename `Sample` -> `biospecimen_id` and `FusionName` -> `final_alteration`
  #  to enable later data joins
  dplyr::rename(biospecimen_id = Sample, 
                final_alteration = FusionName) 



# now can remove the reciprocal: 
rm(fusion_reciprocal_df)

## Prepare the CNV data --------------------------------------

cnv_df <- cnv_df %>%
  # filter to only relevant genes and ids
  dplyr::filter(gene_symbol %in% gene_list,
                biospecimen_id %in% bs_ids) %>%
  # create a single variable indicating the type of CNV we have
  dplyr::mutate(
    final_alteration = glue::glue("{cytoband}-{status}-{copy_number}")
  ) %>%
  # select only relevant columns moving forward
  dplyr::select(
    biospecimen_id, 
    Hugo_Symbol = gene_symbol, 
    final_alteration
  )



## Combine alterations into single data frame and export -------------------------------

# At this point, the three data frames have the same three columns:
# biospecimen_ID, Hugo_Symbol, final_alteration

alteration_df <- dplyr::bind_rows(
  maf_df, 
  fusion_df, 
  cnv_df) %>%
  # ensure we have rows for all combinations of samples and genes
  tidyr::complete(biospecimen_id, Hugo_Symbol) 


# Now, we need handle situations with >1 alteration for a given gene. 
# We'll group as semi-colon separated list
alteration_df %>%
  dplyr::group_by(biospecimen_id, Hugo_Symbol) %>%
  # define the final.final alteration for this sample at this gene
  dplyr::summarize(final_final_alteration = paste(final_alteration, collapse = "; ")) %>%
  dplyr::ungroup() %>%
  # spread: Genes should be columns, and alterations should be values with `NA` for nothing detected
  # first, group on biospecimen_id so we have _one row_ per biospecimen_id
  dplyr::group_by(biospecimen_id) %>%
  tidyr::spread(Hugo_Symbol, final_final_alteration) %>%
  dplyr::ungroup() %>%
  # rename for joining with metadata
  dplyr::rename(Kids_First_Biospecimen_ID = biospecimen_id) %>%
  # bring in some metadata, using `full_join` to ensure all ids make it in the end
  dplyr::full_join(histologies_df) %>%
  # just in case, though there do not appear to be any duplicate rows!
  dplyr::distinct() %>%
  # rearrange columns a bit
  dplyr::select(
    sample_id, 
    Kids_First_Biospecimen_ID, 
    experimental_strategy, 
    germline_sex_estimate,
    dplyr::contains("display"),
    dplyr::everything()
  ) %>%
  # Finally, arrange on sample_id
  dplyr::arrange(sample_id) 


# Check final.final count:
if (length(unique(zenodo_df$sample_id)) != openpbta_total_samples) {
  stop("An error occurred in final steps.")
} 


# Export to CSV file :tada:
readr::write_csv(zenodo_df, zenodo_csv_file)
  
  
  