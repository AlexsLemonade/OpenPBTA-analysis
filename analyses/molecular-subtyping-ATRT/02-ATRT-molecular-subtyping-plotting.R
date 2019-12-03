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

# Install and load in ComplexHeatMap
if (!("ComplexHeatmap" %in% installed.packages())) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)

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

# Read in RNA expression data
zscored_expression <-
  readr::read_rds(file.path(
    results_dir,
    "atrt_zscored_expression.RDS"
  ))

# Read in final output data.frame from `01-ATRT-molecular-subtyping-data-prep.Rmd`
final_df <-
  readr::read_tsv(file.path(results_dir, "ATRT_molecular_subtypes.tsv.gz"))

#### Plot initial heatmap ----------------------------------

# Read in initial heatmap data.frame
initial_ha_atrt_df <-
  readr::read_rds(file.path(results_dir, "initial_heatmap_annotation.RDS"))

# Create a correlation matrix of the expression data
initial_mat <- cor(zscored_expression, method = "pearson")

# Make the initial annotation data.frame a Heatmap Annotation object
initial_ha_atrt <- HeatmapAnnotation(df = initial_ha_atrt_df)

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

# Create a numeric expression matrix for final heatmap
atrt_expression_mat <- final_df %>%
  dplyr::select(
    -c(
      "sample_id",
      "Kids_First_Participant_ID",
      "location_summary",
      "primary_site",
      "age_at_diagnosis_days",
      "germline_sex_estimate",
      "SMARCB1_focal_status",
      "SMARCA4_focal_status",
      "tmb"
    )
  ) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  t()

# Create a correlation matrix
final_mat <- cor(atrt_expression_mat, method = "pearson")

# Read in final heatmap annotation data.frame
final_ha_atrt_df <-
  readr::read_rds(file.path(results_dir, "final_heatmap_annotation.RDS"))

# Make the final annotation data.frame a Heatmap Annotation object
final_ha_atrt <- HeatmapAnnotation(df = final_ha_atrt_df)

# Plot and save heatmap
png(
  file.path(plots_dir, "final_annotated_heatmap.png"),
  width = 820,
  height = 604,
  units = "px"
)
Heatmap(
  final_mat,
  heatmap_legend_param = list(title = "correlation"),
  top_annotation = final_ha_atrt
)
dev.off()
