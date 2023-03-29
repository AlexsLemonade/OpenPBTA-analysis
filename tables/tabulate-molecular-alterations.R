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

  
# Source util functions for data processing:
source(
  file.path(root_dir, "tables", "util", "functions-molecular-alterations.R")
)


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



#### Prepare maf, cnv, fusion data ------------------

maf_df <- prepare_maf(maf_df)
fusion_df <- prepare_fusion(fusion_df)
cnv_df <- prepare_cnv(cnv_df)


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
if (!(length(unique(alterations_df$Hugo_Symbol)) == length(gene_list))) {
  stop("Some genes are missing from the alteration df.")
}

# Convert to presence/absence of _all combinations_ of genes and variants for each sample_id
alterations_presence_df <- alterations_df %>%
  # If the row is here already, then the alteration is present
  dplyr::mutate(present = "Y") %>% 
  # Fill in all combinatins of gene/variant for each sample, 
  #  filling in "N" (no) for present for any new rows 
  tidyr::complete(sample_id, Hugo_Symbol, Variant_Classification, 
                  fill = list(present = "N")) %>%
  dplyr::distinct()

  
# At this point, alterations_presence_df is a Y/N tibble, but it contains only `sample_id`, not DNA and RNA biospecimen ids.
# The final step is to add that information in, but we need to handle different groups separately:
# - First, the fully non-ambiguous samples
# - Second, the ambiguous samples: where DNA and RNA ids cannot be reliably linked, OR 
#   where there are multiple same-modality experiments (since we are working with all samples, not `primary-only`)

## Combine alterations_presence_df with biospecimen sample ids ----------------------
  
## First the fully non-ambiguous samples:
export_identifiable_ids <- tumor_sample_ids_df %>%
  # only ids of interest
  dplyr::filter(!(ambiguous), 
                !(multiple_exp_strategy)) %>%
  # spread ids to get separate RNA and DNA biospecimen columns:
  tidyr::spread(modality, Kids_First_Biospecimen_ID) %>%
  # join with alterations themselves
  dplyr::inner_join(
    alterations_presence_df,
    by = "sample_id"
  ) %>%
  # update columns
  dplyr::select(
    sample_id, 
    Kids_First_Biospecimen_ID_DNA = DNA,
    Kids_First_Biospecimen_ID_RNA = RNA,
    broad_histology_display,
    cancer_group_display,
    germline_sex_estimate,
    Hugo_Symbol,
    Variant_Classification,
    present
  ) %>%
  # finally, enwiden:
  tidyr::spread(Variant_Classification, present)


# Next, the ambiguous ids:
#........


  
  
  
  