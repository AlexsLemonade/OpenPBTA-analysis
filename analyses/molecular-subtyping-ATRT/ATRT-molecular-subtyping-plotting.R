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
# Rscript 'analyses/molecular-subtyping-ATRT/ATRT-molecular-subtyping-plotting.R'

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# Load in ComplexHeatMap
if (!("ComplexHeatmap" %in% installed.packages())) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)
if (!("circlize" %in% installed.packages())) {
  install.packages("circlize")
}

# Set seed for Heatmap function
set.seed(2019)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results and plots directory
results_dir <-
  file.path(root_dir, "scratch")

plots_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Source script with data preparation
source(
  file.path(
    root_dir,
    "analyses",
    "molecular-subtyping-ATRT",
    "ATRT-molecular-subtyping.R"
  )
)

#### Custom function -----------------------------------------------------------

get_mode <- function(x){
  return(names(sort(table(x), decreasing = T, na.last = T)[1]))
}

#### Plot initial heatmap ------------------------------------------------------

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
    Kids_First_Participant_ID,
    primary_site,
    age_at_diagnosis_days
  ))

# Make the annotation data.frame a Heatmap Annotation object
intial_ha_atrt <- HeatmapAnnotation(df = initial_annotation_df)

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
  dplyr::select(
    Kids_First_Biospecimen_ID,
    z_score,
    MYC_V1_ssgsea_score,
    MYC_V2_ssgsea_score,
    Notch_ssgsea_score
  ) %>%
  dplyr::group_by(Kids_First_Biospecimen_ID) %>%
  dplyr::summarise(
    z_score = mean(!is.na(z_score)),
    MYC_V1_ssgsea_score = mean(MYC_V1_ssgsea_score),
    MYC_V2_ssgsea_score = mean(MYC_V2_ssgsea_score),
    Notch_ssgsea_score = mean(Notch_ssgsea_score)
  ) %>%
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
  dplyr::group_by(Kids_First_Biospecimen_ID) %>%
  dplyr::distinct() %>%
  dplyr::summarise(
    SMARCB1_focal_status = get_mode(SMARCB1_focal_status),
    location_summary = get_mode(location_summary),
    SMARCA4_focal_status = get_mode(SMARCA4_focal_status),
    Overexpressed_gene_sets = get_mode(Overexpressed_gene_sets)
    # Used the `get_mode` custom function above to solve the
    # `Column must be length 1 (a summary value)` error, but not entirely sure
    # this is the right fix here
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
