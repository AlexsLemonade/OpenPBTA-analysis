# This script subsets the focal copy number, RNA expression, tumor mutation 
# burden and histologies` files to include only ATRT samples.

# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/00-subset-files-for-ATRT.R'


#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results and plots directory
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "atrt-subset")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

# Select wanted columns in metadata for merging and assign to a new object
select_metadata <- metadata %>%
  dplyr::select(sample_id,
                Kids_First_Participant_ID,
                Kids_First_Biospecimen_ID,
                glioma_brain_region)

# Read in ssGSEA pathway information
ssGSEA <-
  as.data.frame(readr::read_rds(
    file.path(
      root_dir,
      "analyses",
      "ssgsea-hallmark",
      "results",
      "GeneSetExpressionMatrix.RDS"
    )
  ))

# Read in RNA expression data
stranded_expression <-
  readr::read_rds(
    file.path(
      root_dir,
      "data",
      "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
    )
  )

# Read in focal CN data
## TODO: This section will be updated to read in focal CN data derived from
##       copy number consensus calls.
cn_df <- readr::read_tsv(
  file.path(
    root_dir,
    "analyses",
    "focal-cn-file-preparation",
    "results",
    "controlfreec_annotated_cn_autosomes.tsv.gz"
  )
)

# Read in consensus mutation data
tmb_df <-
  data.table::fread(file.path(root_dir,
                              "data",
                              "pbta-snv-consensus-mutation-tmb.tsv"))

#### Filter metadata -----------------------------------------------------------

# Define regions of the brain (using Anatomy of the Brain figure found at
# https://www.ncbi.nlm.nih.gov/books/NBK65903/figure/CDR0000574573__205/)
supratentorial <-
  c(
    "Skull",
    "Thalamus",
    "Temporal Lobe",
    "Frontal Lobe",
    "Parietal Lobe",
    "Cerebrum",
    "Basal Ganglia",
    "Cranial Nerves NOS",
    "Basal Ganglia;Temporal Lobe",
    "Frontal Lobe;Parietal Lobe;Temporal Lobe",
    "Parietal Lobe;Temporal Lobe",
    "Frontal Lobe;Parietal Lobe"
  )

infratentorial <-
  c(
    "Cerebellum/Posterior Fossa",
    "Brain Stem- Pons;Cerebellum/Posterior Fossa",
    "Cerebellum/Posterior Fossa;Other locations NOS",
    "Brain Stem",
    "Brain Stem- Midbrain/Tectum;Ventricles",
    "Cerebellum/Posterior Fossa;Ventricles",
    "Cerebellum/Posterior Fossa;Spinal Cord- Cervical;Spinal Cord- Lumbar/Thecal Sac;Spinal Cord- Thoracic",
    "Other locations NOS;Spinal Cord- Lumbar/Thecal Sac;Spinal Cord- Thoracic;Ventricles"
  )

# Filter metadata for `ATRT` and define `location_summary` based on values in
# `primary_site`
atrt_df <- metadata %>%
  dplyr::filter(short_histology == "ATRT",
                experimental_strategy == "RNA-Seq") %>%
  dplyr::mutate(
    location_summary = dplyr::case_when(
      primary_site %in% infratentorial ~ "infratentorial",
      primary_site %in% supratentorial ~ "supratentorial",
      TRUE ~ "NA"
    )
  ) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarize(
    Kids_First_Biospecimen_ID = paste(sort(unique(
      Kids_First_Biospecimen_ID
    )),
    collapse = ", "),
    Kids_First_Participant_ID,
    location_summary,
    age_at_diagnosis_days,
    germline_sex_estimate,
    primary_site
  )

# Write to file
readr::write_tsv(atrt_df, file.path(results_dir, "atrt_histologies.tsv"))

#### Filter expression data ----------------------------------------------------

# Filter to ATRT samples only -- we can use atrt_df because it is subset to
# RNA-seq samples
stranded_expression <- stranded_expression %>%
  dplyr::select(atrt_df$Kids_First_Biospecimen_ID)

# Log2 transformation
norm_expression <- log2(stranded_expression + 1)

# normData mean and sd
norm_expression_means <- rowMeans(norm_expression, na.rm = TRUE)
norm_expression_sd <- apply(norm_expression, 1, sd, na.rm = TRUE)

# Subtract mean
expression_zscored <-
  sweep(norm_expression, 1, norm_expression_means, FUN = "-")

# Divide by SD remove NAs and Inf values from zscore for genes with 0 in normData
expression_zscored <-
  sweep(expression_zscored, 1, norm_expression_sd, FUN = "/") %>%
  dplyr::na_if(Inf) %>%
  na.omit()

# Save matrix with all genes to file for downstream plotting
readr::write_rds(expression_zscored, 
                 file.path(results_dir, "atrt_zscored_expression.RDS"))

#### Filter focal CN data ------------------------------------------------------

# Filter focal CN to ATRT samples only 
cn_metadata <- cn_df %>%
  dplyr::left_join(select_metadata,
                   by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(gene_symbol,
                sample_id,
                Kids_First_Participant_ID,
                biospecimen_id,
                status) %>%
  dplyr::distinct()

# Write to file
readr::write_tsv(cn_metadata, file.path(results_dir, "atrt_focal_cn.tsv.gz"))

#### Filter ssGSEA data --------------------------------------------------------

# Calculate ssGSEA mean and sd
ssGSEA_means <- rowMeans(ssGSEA, na.rm = TRUE)
ssGSEA_sd <- apply(ssGSEA, 1, sd, na.rm = TRUE)

# Subtract mean
ssGSEA_zscored <- sweep(ssGSEA, 1, ssGSEA_means, FUN = "-")

# Divide by SD remove NAs and Inf values from zscore for genes with 0
ssGSEA_zscored <-
  sweep(ssGSEA_zscored, 1, ssGSEA_sd, FUN = "/") %>%
  dplyr::na_if(Inf) %>%
  na.omit()

# Transpose
transposed_ssGSEA <- t(ssGSEA_zscored)

# Select wanted pathways and merge metadata
transposed_ssGSEA <- transposed_ssGSEA %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    HALLMARK_MYC_TARGETS_V1 = mean(HALLMARK_MYC_TARGETS_V1),
    HALLMARK_MYC_TARGETS_V2 = mean(HALLMARK_MYC_TARGETS_V2),
    HALLMARK_NOTCH_SIGNALING = mean(HALLMARK_NOTCH_SIGNALING)
  )

# Write to file
readr::write_tsv(transposed_ssGSEA, file.path(results_dir, "atrt_ssgsea.tsv"))

#### Filter tumor mutation burden data -----------------------------------------

tmb_df <- tmb_df %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::inner_join(select_metadata,
                    by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID"))

# Write to file 
readr::write_tsv(tmb_df, file.path(results_dir, "atrt_tmb.tsv"))
