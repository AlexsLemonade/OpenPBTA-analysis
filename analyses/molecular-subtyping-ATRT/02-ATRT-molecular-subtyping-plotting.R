# This script addresses the issue of molecular subtyping ATRT samples by 
# plotting the filtered ATRT data produced in the `ATRT-molecular-subtyping.R`
# script.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/02-ATRT-molecular-subtyping-plotting.R'

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# Load in ComplexHeatMap
if (!("ComplexHeatmap" %in% installed.packages())) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)

# Set seed for Heatmap function
set.seed(2019)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results and plots directory
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "results")

plots_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

# Read in RNA expression data
stranded_expression <-
  readr::read_rds(file.path(
    root_dir,
    "data",
    "pbta-gene-expression-rsem-fpkm.stranded.rds"
  ))

# Read in output file produced in `01-ATRT-molecular-subtyping-data-prep.Rmd`
atrt_expression_cn_tmb_df <-
  readr::read_tsv(file.path(results_dir, "ATRT_molecular_subtypes.tsv.gz"))

#### Filter for ATRT and plot initial heatmap ----------------------------------

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

# Filter metadata for ATRT samples
atrt_df <- metadata %>%
  dplyr::filter(short_histology == "ATRT",
                experimental_strategy == "RNA-Seq") %>%
  dplyr::mutate(
    location_summary = dplyr::case_when(
      primary_site %in% infratentorial ~ "infratentorial",
      primary_site %in% supratentorial ~ "supratentorial",
      TRUE ~ as.character(primary_site)
    )
  ) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarize(Kids_First_Biospecimen_ID = paste(sort(unique(
    Kids_First_Biospecimen_ID
  )),
  collapse = ", "),
  Kids_First_Participant_ID,
  location_summary,
  age_at_diagnosis_days,
  germline_sex_estimate, 
  primary_site) 

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

# Create a correlation matrix of the expression data
initial_mat <- cor(stranded_expression_df, method = "pearson")

# Create an annotation data.frame for the relevant annotation data
initial_annotation_df <- atrt_df %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  dplyr::select(-c(
    sample_id,
    Kids_First_Participant_ID,
    primary_site,
    age_at_diagnosis_days
  ))

# Make the annotation data.frame a Heatmap Annotation object
initial_ha_atrt <- HeatmapAnnotation(df = initial_annotation_df)

# Plot and save the Heatmap
png(
  file.path(plots_dir, "initial_heatmap.png"),
  width = 820,
  height = 604,
  units = "px"
)
Heatmap(
  initial_mat,
  heatmap_legend_param = list(title = "correlation"),
  top_annotation = initial_ha_atrt
)
dev.off()

#### Plot final annotated heatmap ----------------------------------------------

# Create a numeric expression matrix for plotting
atrt_expression_mat <- atrt_expression_cn_tmb_df %>%
  dplyr::select(-c("sample_id",
                   "Kids_First_Participant_ID",
                   "location_summary",
                   "age_at_diagnosis_days",
                   "germline_sex_estimate",
                   "SMARCB1_focal_status",
                   "SMARCA4_focal_status",
                   "Overexpressed_gene_sets",
                   "tmb"
  )) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  t()

atrt_expression_mat <- cor(atrt_expression_mat, method = "pearson")

# Create an annotation matrix for the atrt expression data
annotation_mat <- atrt_expression_cn_tmb_df %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    location_summary,
    SMARCB1_focal_status,
    SMARCA4_focal_status,
    Overexpressed_gene_sets
  ) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# Make the annotation matrix a Heatmap Annotation object
ha_atrt <- HeatmapAnnotation(df = annotation_mat)

# Plot and save heatmap
png(
  file.path(plots_dir, "final_annotated_heatmap.png"),
  width = 820,
  height = 604,
  units = "px"
)
Heatmap(
  atrt_expression_mat,
  heatmap_legend_param = list(title = "correlation"),
  top_annotation = ha_atrt
)
dev.off()
