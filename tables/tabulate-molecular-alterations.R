# S. Spielman 2023 for CCDL
# 
# This script tabulates, for _all samples_ in OpenPBTA cohort, the presence/absence of 
#  alterations depicted in the manuscript's oncoprint plots.
# This script is heavily inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R
# 
# Output file structure:
#  sample_id    DNA BS     RNA BS     mapping_type       <------ alterations ------>
#  XXXX-YYYY    BS_###     BS_###    <1:1 or 1:many>            Y/N values
#
# Note that the oncoprint module is assumed to have been run before this script, 
#  such that `scratch/oncoprint-plots/` is populated.

# Define pipe
`%>%` <- dplyr::`%>%`

# helper variable for later - there should be 1074 tumors
openpbta_n_tumors <- 1074

# Set up paths ------------
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

### Helper functions ----------------------------
# This section holds functions for preparing maf, cnv, and fusion data so they can be
#  combined into a single data frame containing all alterations in the given tumors.
# These functions are **NOT TO BE USED** outside of this script/context.
# This code is heavily/entirely inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R



prepare_maf <- function(maf_df, tumor_sample_ids_df) {
  maf_df %>%
    dplyr::inner_join(
      tumor_sample_ids_df,
      by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")
    ) %>%
    # now let's remove this `Tumor_Sample_Barcode` column with biospecimen IDs in
    # preparation for our next step -- renaming `sample_id`
    dplyr::select(-Tumor_Sample_Barcode) %>%
    dplyr::rename(Tumor_Sample_Barcode = sample_id)
}



prepare_fusion <- function(fusion_df, tumor_sample_ids_df) {
  
  # We'll handle fusions where reciprocal fusions exist (e.g., 
  # reciprocal_exists == TRUE) separately from other fusions
  fusion_reciprocal_df <- fusion_df %>%
    # Reciprocal fusions only
    dplyr::filter(reciprocal_exists) %>%
    # BSID + Gene1--Gene2
    dplyr::select(Sample, FusionName) %>%
    # Because we're only looking at presence or absence here, we can filter to 
    # distinct identifier-fusion pairs
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
    # When fusions are in alphabetical order, this ensures each reciprocal
    # fusion is only counted once
    dplyr::distinct()
  
  # No need to put fusions where no reciprocal exists in alphabetical order,
  # so we handle these separately
  fusion_no_reciprocal_df <- fusion_df %>%
    dplyr::filter(!reciprocal_exists) %>%
    dplyr::select(Sample, FusionName) %>%
    # But we can remove duplicates to avoid them being counted as multihit
    dplyr::distinct()
  
  # Bind reciprocal and no reciprocal together
  fusion_filtered_df <- dplyr::bind_rows(fusion_reciprocal_df, fusion_no_reciprocal_df)
  rm(fusion_reciprocal_df, fusion_no_reciprocal_df)
  
  # Separate fusion gene partners
  fus_sep <- fusion_filtered_df %>%
    # Separate the 5' and 3' genes
    tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--") %>%
    # Use row numbers to mark unique fusions - this will help us when
    # we melt and remove selfie fusions below
    tibble::rowid_to_column("Fusion_ID") %>%
    dplyr::select(Fusion_ID, Sample, Gene1, Gene2) %>%
    reshape2::melt(id.vars = c("Fusion_ID", "Sample"),
                   variable.name = "Partner",
                   value.name = "Hugo_Symbol") %>%
    dplyr::arrange(Fusion_ID)
  
  multihit_fusions <- fus_sep %>%
    # For looking at multi-hit fusions, we want to remove selfie fusions
    # If we drop the 5', 3' information in Partner, we can use distinct
    # to remove the selfie fusions because Fusion_ID marks what fusion
    # it came from
    dplyr::select(-Partner) %>%
    dplyr::distinct() %>%
    # Now we want to count how many times a gene is fused within a sample
    # Anything with more than one count is considered multi-hit
    dplyr::count(Sample, Hugo_Symbol) %>%
    dplyr::filter(n > 1) %>%
    dplyr::select(-n) %>%
    dplyr::mutate(Variant_Classification = "Multi_Hit_Fusion")
  
  # Filter out multi-hit fusions from the other fusions that we will
  # label based on whether they are the 5' or 3' gene
  single_fusion <- fus_sep %>%
    dplyr::anti_join(multihit_fusions,
                     by = c("Sample", "Hugo_Symbol")) %>%
    dplyr::select(Sample, Hugo_Symbol) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Variant_Classification = "Fusion")
  
  # Combine multi-hit and single fusions into one data frame
  # And add the other identifiers!
  reformat_fusion_df <- multihit_fusions %>%
    dplyr::bind_rows(single_fusion) %>%
    dplyr::mutate(Variant_Type = "OTHER") %>%
    dplyr::inner_join(tumor_sample_ids_df,
                      by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
    dplyr::rename(Tumor_Sample_Barcode = sample_id,
                  Kids_First_Biospecimen_ID = Sample)
  
  # return reformatted
  reformat_fusion_df
}



prepare_cnv <- function(cnv_df, tumor_sample_ids_df) {
  
  cnv_df %>%
    dplyr::inner_join(
      tumor_sample_ids_df,
      by = c("biospecimen_id"="Kids_First_Biospecimen_ID")) %>%
    dplyr::mutate(Tumor_Sample_Barcode =  sample_id) %>%
    dplyr::rename(Variant_Classification = status,
                  Hugo_Symbol = gene_symbol) %>%
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
    # mutate loss and amplification to Del and Amp to fit Maftools format
    dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "deep deletion" ~ "Del",
                                                            Variant_Classification == "amplification" ~ "Amp",
                                                            TRUE ~ as.character(Variant_Classification))) %>%
    # only keep Del and Amp calls
    dplyr::filter(Variant_Classification %in% c("Del", "Amp"))
  
}


#### Read in data ------------------------------------
histologies_df <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  dplyr::inner_join(
    readr::read_tsv(pal_file) %>%
      dplyr::select(cancer_group, cancer_group_display, broad_histology, broad_histology_display)
  )

# Only keep required maf columns 
keep_cols <-c("Hugo_Symbol", 
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2",
              "Variant_Classification",
              "Variant_Type",
              "Tumor_Sample_Barcode",
              "HGVSp_Short")

maf_df <- readr::read_tsv(maf_file) %>%
  dplyr::select(keep_cols)

hotspots_maf_df <- readr::read_tsv(hotspots_maf_file) %>%
  dplyr::select(keep_cols)
  
# merge hotspots maf to input maf
maf_df <- maf_df %>%
  dplyr::bind_rows(hotspots_maf_df) %>%
  dplyr::distinct()


# read fusion data
fusion_df <- readr::read_tsv(fusion_file)


# read and join cnv with histology information
cnv_autosomes_df <- readr::read_tsv(cnv_autosomes_file) %>%
  dplyr::left_join(
    dplyr::select(histologies_df, 
                  Kids_First_Biospecimen_ID,
                  germline_sex_estimate),
            by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")
  )

# combine autosomes with x/y cnv
cnv_df <- dplyr::bind_rows(
  cnv_autosomes_df, 
  readr::read_tsv(cnv_xy_file)
)


## Filter to relevant tumors --------------------------------------

# The final output should only tabuluate alterations for N=1074 tumors 
# We'll need to identify those samples, including their sample_id and biospecimen.

# In addition, we'll note if a sample_id is ambiguous. Ambiguous sample_id's will 
# have more than 2 rows associated with it in the histologies file when looking at 
# tumor samples -- that means we won't necessarily be able to determine **when an 
# WGS/WXS assay maps to an RNA-seq assay** for the purpose of tabulating molecular alterations

tumor_sample_ids_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue") %>%
  # keep these columns of interest moving forward:
  dplyr::select(sample_id, 
                Kids_First_Biospecimen_ID, 
                cancer_group_display, 
                broad_histology_display, 
                experimental_strategy, 
                germline_sex_estimate) %>%
  dplyr::distinct() %>%
  # annotate if it is ambiguous or: any rows with >2 sample_ids
  # Associated with 2 sample_ids that have the same experimental strategy
  dplyr::add_count(sample_id) %>%
  dplyr::mutate(ambiguous = n > 2) %>% 
  # now, annotate if it has multiple rows of SAME experimental_strategy (needed for later wrangling)
  dplyr::add_count(sample_id, experimental_strategy) %>%
  dplyr::mutate(multiple_exp_strategy = n >= 2) %>% 
  # remove temporary counting column
  dplyr::select(-n) %>%
  # add a modality column for later bookkeeping
  dplyr::mutate(modality = dplyr::case_when(
    experimental_strategy == "RNA-Seq" ~ "RNA",
    experimental_strategy == "WGS"     ~ "DNA",
    experimental_strategy == "WXS"     ~ "DNA",
    experimental_strategy == "Targeted Sequencing" ~ "DNA"
  )) %>%
  # remove experimental_strategy to enable later joining
  dplyr::select(-experimental_strategy)


# Check that this is equal to `openpbta_n_tumors`
if (!(length(unique(tumor_sample_ids_df$sample_id))) == openpbta_n_tumors) {
  stop("Bad sample_id parsing.")
} 

# Pull out biospecimen id column for convenience 
bs_ids <- tumor_sample_ids_df %>%
  dplyr::pull(Kids_First_Biospecimen_ID) %>%
  unique() 

# Filter dfs to those ids:
histologies_df <-  dplyr::filter(histologies_df, Kids_First_Biospecimen_ID %in% bs_ids) 
maf_df         <-  dplyr::filter(maf_df, Tumor_Sample_Barcode %in% bs_ids) 
fusion_df      <-  dplyr::filter(fusion_df, Sample %in% bs_ids) 
cnv_df         <-  dplyr::filter(cnv_df, biospecimen_id %in% bs_ids) 

# What are all the bs ids that end up being considered? 
# This variable can be used when joining things back up
final_bs_ids <- length(unique(
  c(
    maf_df$Tumor_Sample_Barcode,
    fusion_df$Sample,
    cnv_df$biospecimen_id
)))



#### Prepare maf, cnv, fusion data to be combined ------------------

maf_df <- prepare_maf(maf_df, tumor_sample_ids_df)
fusion_df <- prepare_fusion(fusion_df, tumor_sample_ids_df)
cnv_df <- prepare_cnv(cnv_df, tumor_sample_ids_df)


### Combine alterations and filter to GOI ---------------------

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

select_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")


# Combine all alterations:
alterations_df <- list(maf_df, fusion_df, cnv_df) %>%
  # keeping it interesting
  purrr::map(
    dplyr::select,
    select_cols
  ) %>%
  dplyr::bind_rows() %>%
  # rename Tumor_Sample_Barcode to sample_id
  dplyr::rename(sample_id = Tumor_Sample_Barcode) %>%
  # filter to only genes in gene_list
  dplyr::filter(Hugo_Symbol %in% gene_list) %>%
  # some wrangling to help with dev:
  dplyr::arrange(Hugo_Symbol) %>%
  tibble::as_tibble()

# Ensure we still have all the genes at least _somewhere_ in the data frame
if (length(unique(alterations_df$Hugo_Symbol)) != length(gene_list)) {
  stop("Some genes are missing from the alteration df.")
}

# Convert to presence/absence of _all combinations_ of genes and variants for each sample_id
alterations_presence_df <- alterations_df %>%
  # If the row is here already, then the alteration is present
  dplyr::mutate(present = TRUE) %>% 
  # Fill in all combinatins of gene/variant for each sample, 
  #  filling in FALSE for present for any new rows 
  tidyr::complete(sample_id, Hugo_Symbol, Variant_Classification, 
                  fill = list(present = FALSE)) %>%
  dplyr::distinct()

  
# At this point, alterations_presence_df is a Y/N tibble, but it contains only `sample_id`, not DNA and RNA biospecimen ids.
# The final step is to add that information in, but we need to handle different groups separately:
# - First, the fully non-ambiguous samples
# - Second, the ambiguous samples: where DNA and RNA ids cannot be reliably linked, OR 
#   where there are multiple same-modality experiments (since we are working with all samples, not `primary-only`)

## Combine alterations_presence_df with biospecimen sample ids ----------------------
  
# First, let's "combine" the ambiguous ids: 
# For multiple DNA or RNA for a given sample_id, convert to a semi-colon separated single value
# For example, two DNA ids on separate rows named BS_1 and BS_2 would be collapsed into `BS_1; BS2`
merged_ambiguous_ids <- tumor_sample_ids_df %>%
  # get only samples of interest
  dplyr::filter(ambiguous | multiple_exp_strategy) %>%
  # for each sample_id, combine biospecimen ids by modality to have 
  # at most one "value" (semi-colon separated list) for each modality 
  dplyr::group_by(sample_id, modality) %>%
  dplyr::summarize(Kids_First_Biospecimen_ID = paste(Kids_First_Biospecimen_ID, collapse="; ")) %>%
  dplyr::ungroup() %>%
  # bring back the other variables from `tumr_sample_ids_df`
  dplyr::inner_join(
    dplyr::select(
      tumor_sample_ids_df, 
      sample_id, 
      broad_histology_display,
      cancer_group_display,
      germline_sex_estimate
    )
  ) %>%
  dplyr::distinct() %>%
  # add indicator that these are ambiguous mappings
  dplyr::mutate(ambiguous_biospecimen_mapping = TRUE)


# Now, subset to just the identifiable ids and bind_rows with `merged_ambiguous_ids` for final processing
prepared_ids_df <- tumor_sample_ids_df %>%
  # only ids of interest
  dplyr::filter(!(ambiguous), 
                !(multiple_exp_strategy)) %>%
  dplyr::select(-ambiguous, -multiple_exp_strategy) %>%
  # add indicator that these are NOT ambiguous mappings
  dplyr::mutate(ambiguous_biospecimen_mapping = FALSE) %>%
  # combine with prepared ambiguous samples
  dplyr::bind_rows(merged_ambiguous_ids) %>%
  dplyr::distinct()

# Again, check that we've got everything: there should be `openpbta_n_tumors` 
#  unique sample_id values
if (!(length(unique(tumor_sample_ids_df$sample_id))) == openpbta_n_tumors) {
  stop("An error occurred while dealing with ambiguous vs. identifiable biospecimen ids.")
} 

## Finally, process the full data for export:
zenodo_df <- prepared_ids_df %>%
  tidyr::spread(modality, Kids_First_Biospecimen_ID) %>%
  # join with alterations themselves
  dplyr::inner_join(
    alterations_presence_df,
    by = "sample_id"
  ) %>%
  # update columns and their order
  dplyr::select(
    sample_id, 
    Kids_First_Biospecimen_ID_DNA = DNA,
    Kids_First_Biospecimen_ID_RNA = RNA,
    ambiguous_biospecimen_mapping,
    broad_histology_display,
    cancer_group_display,
    germline_sex_estimate,
    Hugo_Symbol,
    Variant_Classification,
    present
  ) %>%
  # finally, enwiden:
  tidyr::spread(Variant_Classification, present)

  
# Export to CSV file :tada:
readr::write_csv(zenodo_df, zenodo_csv_file)
  
  
  