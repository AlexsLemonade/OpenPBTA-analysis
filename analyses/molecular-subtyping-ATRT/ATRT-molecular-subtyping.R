# This script addresses the issue of molecular subtyping ATRT samples.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/ATRT-molecular-subtyping.R'

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# File path to plots directory
plots_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# File path to results directory
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

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
    "Frontal Lobe;Parietal Lobe",
    "Frontal Lobe;Meninges/Dura;Parietal Lobe;Spinal Cord- Cervical;Spinal Cord- Lumbar/Thecal Sac;Spinal Cord- Thoracic;Temporal Lobe;Ventricles"
    # Not sure how the line above should be broken up as `Frontal Lobe` is
    # supratentorial, while `Spinal Cord` is not
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
    "Other locations NOS;Spinal Cord- Lumbar/Thecal Sac;Spinal Cord- Thoracic;Ventricles",
    "Other locations NOS;Spinal Cord- Lumbar/Thecal Sac;Spinal Cord- Thoracic;Temporal Lobe"
    # Not sure how the line above should be broken up as `Temporal Lobe` is
    # supratentorial, while the `Spinal Cord` is not
  )

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

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
  readr::read_rds(file.path(
    root_dir,
    "data",
    "pbta-gene-expression-rsem-fpkm.stranded.rds"
  ))

# Read in focal CN data
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
  data.table::fread(
    file.path(
      root_dir,
      "data",
      "snv-consensus_11122019",
      "consensus_mutation_tmb.tsv"
    )
  )

#### Filter data and plot initial heatmap --------------------------------------

# Filter metadata for `ATRT` and define `location_summary` based on values in
# `primary_site`
atrt_df <- metadata %>%
  dplyr::filter(short_histology == "ATRT",
                experimental_strategy == "RNA-Seq") %>%
  dplyr::select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    age_at_diagnosis_days,
    germline_sex_estimate,
    primary_site
  ) %>%
  dplyr::mutate(
    location_summary = dplyr::case_when(
      primary_site %in% infratentorial ~ "infratentorial",
      primary_site %in% supratentorial ~ "supratentorial",
      TRUE ~ as.character(primary_site)
    )
  )

# Save Biospecimen Ids in a vector
column_names <- colnames(stranded_expression)

# Filter for Ids relevant to `ATRT`
column_names <- column_names %>%
  as.data.frame() %>%
  dplyr::filter(. %in% atrt_df$Kids_First_Biospecimen_ID)

# Make a data.frame with only the expression values for ATRT samples
stranded_expression_df <- stranded_expression %>%
  as.data.frame() %>%
  dplyr::select(column_names$.)

#### Join data -----------------------------------------------------------------

# Format expression data into long format
long_stranded_expression <- stranded_expression %>%
  tidyr::gather(biospecimen_id, expression_value, -gene_id) %>%
  dplyr::distinct() %>%
  dplyr::mutate(gene_id = gsub(".*\\_", "", gene_id)) # Trim the gene ids to include only the gene symbol

# Calculate log2 expression and z-scores of expression values
long_stranded_expression <- long_stranded_expression %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(log2_exp = log2(expression_value + 1),
                z_score = (log2_exp - mean(log2_exp) / sd(log2_exp)))

expression_metadata <- long_stranded_expression %>%
  dplyr::left_join(metadata, by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(biospecimen_id, Kids_First_Participant_ID, z_score, gene_id)

# Join expression data with metadata filtered for `ATRT`
atrt_expression_df <- atrt_df %>%
  dplyr::inner_join(expression_metadata,
                    by = "Kids_First_Participant_ID") %>%
  dplyr::select(-gene_id)

# Define target overexpressed gene vectors
tyr_genes <-
  paste(c("TYR", "MITF", "DCT", "VEGFA", "DNAH11", "SPEF1"), collapse = "|")
shh_genes <-
  paste(c("MYCN", "GLI2", "CDK6", "ASCL1", "HES5/6", "DLL1/3"),
        collapse = "|")
myc_genes <- paste(c("MYC", "HOTAIR", "HOX"), collapse = "|")

# Join focal CN data with metadata
cn_metadata <- cn_df %>%
  dplyr::left_join(metadata,
                   by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(gene_symbol,
                Kids_First_Participant_ID,
                biospecimen_id,
                status) %>%
  dplyr::filter(Kids_First_Participant_ID %in% atrt_expression_df$Kids_First_Participant_ID) %>%
  dplyr::mutate(
    SMARCB1_focal_status = dplyr::case_when(gene_symbol == "SMARCB1" ~ status,
                                            TRUE ~ "neutral"),
    SMARCA4_focal_status = dplyr::case_when(gene_symbol == "SMARCA4" ~ status,
                                            TRUE ~ "neutral"),
    Overexpressed_gene_sets = dplyr::case_when(
      stringr::str_detect(gene_symbol, tyr_genes) &
        status != "loss" ~ "TYR_genes",
      stringr::str_detect(gene_symbol, shh_genes) &
        status != "loss" ~ "SHH_genes",
      stringr::str_detect(gene_symbol, myc_genes) &
        status != "loss" ~ "MYC_genes",
      TRUE ~ "NA"
    )
  ) %>%
  dplyr::select(-gene_symbol) %>%
  dplyr::distinct()

# Join ATRT expression data with focal CN data
atrt_expression_cn_df <- atrt_expression_df %>%
  dplyr::inner_join(cn_metadata, by = "Kids_First_Participant_ID") %>%
  dplyr::select(-status)

# Transpose and filter ssGSEA pathway data
transposed_ssgsea <- t(ssGSEA) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("biospecimen_id") %>%
  dplyr::filter(biospecimen_id %in% atrt_expression_df$Kids_First_Biospecimen_ID) %>%
  dplyr::select(
    biospecimen_id,
    HALLMARK_MYC_TARGETS_V1,
    HALLMARK_MYC_TARGETS_V2,
    HALLMARK_NOTCH_SIGNALING
  )

# Join ATRT expression and focal CN data with transposed ssGSEA data, calculate
# z-scores for each pathway variable
atrt_expression_cn_df <- atrt_expression_cn_df %>%
  dplyr::inner_join(transposed_ssgsea,
                    by = c("Kids_First_Biospecimen_ID" = "biospecimen_id")) %>%
  dplyr::mutate(
    MYC_V1_ssgsea_score = HALLMARK_MYC_TARGETS_V1 - mean(HALLMARK_MYC_TARGETS_V1) / sd(HALLMARK_MYC_TARGETS_V1),
    MYC_V2_ssgsea_score = HALLMARK_MYC_TARGETS_V2 - mean(HALLMARK_MYC_TARGETS_V2) / sd(HALLMARK_MYC_TARGETS_V2),
    Notch_ssgsea_score = HALLMARK_NOTCH_SIGNALING - mean(HALLMARK_NOTCH_SIGNALING) / sd(HALLMARK_NOTCH_SIGNALING)
  ) %>%
  dplyr::select(
    -c(
      "HALLMARK_MYC_TARGETS_V1",
      "HALLMARK_MYC_TARGETS_V2",
      "HALLMARK_NOTCH_SIGNALING"
    )
  )

# Join tumor mutuation data with metadata
tmb_df <- tmb_df %>%
  dplyr::left_join(metadata,
                   by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Tumor_Sample_Barcode, Kids_First_Participant_ID, tmb)

atrt_expression_cn_tmb_df <- atrt_expression_cn_df %>%
  dplyr::left_join(tmb_df, by = "Kids_First_Participant_ID")

## TODO: Add a column to this data.frame denoting `chr22q` loss using the SV
# data.

#### Save results --------------------------------------------------------------

# Save final data.frame
final_df <- atrt_expression_cn_tmb_df %>%
  dplyr::select(-c(
    "biospecimen_id.x",
    "biospecimen_id.y",
    "Tumor_Sample_Barcode"
  )) %>%
  dplyr::distinct() %>%
  dplyr::arrange(desc(tmb))

readr::write_tsv(final_df,
                 file.path(results_dir, "ATRT_molecular_subtypes.tsv.gz"))
