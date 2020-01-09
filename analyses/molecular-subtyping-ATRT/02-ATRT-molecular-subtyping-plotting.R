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

# Install and load in ComplexHeatMap
if (!("ComplexHeatmap" %in% installed.packages())) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)

if (!("matrixStats" %in% installed.packages())) {
  install.packages("matrixStats")
}

if (!("ggfortify" %in% installed.packages())) {
  install.packages("ggfortify")
}

library(ggfortify)
library(dplyr)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results, plots, and subset files directories
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-ATRT")
results_dir <- file.path(module_dir, "results")
plots_dir <- file.path(module_dir, "plots")
subset_dir <- file.path(module_dir, "atrt-subset")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in RNA expression data
log_expression <-
  readr::read_rds(file.path(
    subset_dir,
    "atrt_log_expression.RDS"
  ))

# Read in the subset histologies file
metadata_df <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# Read in final output data.frame from `01-ATRT-molecular-subtyping-data-prep.Rmd`
final_df <-
  readr::read_tsv(file.path(results_dir, "ATRT_molecular_subtypes.tsv"))

#### Plot heatmap --------------------------------------------------------------

# Find and filter to high variance genes
row_variances <- matrixStats::rowVars(as.matrix(log_expression))
high_var_exp <- log_expression[which(row_variances > quantile(row_variances, 0.8)), ]

# Z-score the genes and change the identifiers
high_var_exp <- scale(t(high_var_exp),  # scale works on columns
                      center = TRUE,
                      scale = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  # Convert from the biospecimen IDs to the sample_id -- this will allow us
  # to annotate the heatmap we will make
  inner_join(
    select(
      metadata_df,
      Kids_First_Biospecimen_ID,
      sample_id
    ),
    by = "Kids_First_Biospecimen_ID") %>%
  select(-Kids_First_Biospecimen_ID) %>%
  tibble::column_to_rownames("sample_id") %>%
  # Make it such that rows are genes again
  t()

# Make data.frame that will be used to form an annotation bar for the heatmap
annotation_df <- final_df %>%
  filter(sample_id %in% colnames(high_var_exp)) %>%
  select(
    sample_id,
    location_summary,
    germline_sex_estimate,
    SMARCB1_focal_status,
    SMARCA4_focal_status,
    chr_22q_loss
  ) %>%
  as.data.frame() %>%  # ComplexHeatmap doesn't like tibbles
  tibble::column_to_rownames("sample_id")

# Make into an annotation object
column_annotation <- HeatmapAnnotation(df = annotation_df)

# Plot and save the Heatmap
png(
  file.path(plots_dir, "atrt_heatmap.png"),
  width = 820,
  height = 604,
  units = "px"
)
Heatmap(
  high_var_exp,
  heatmap_legend_param = list(title = "zscore"),
  top_annotation = column_annotation,
  show_row_names = FALSE,
  show_column_names = FALSE
)
dev.off()

# Make a PCA plot
pca_results <- prcomp(t(high_var_exp))
pca_plot <- autoplot(pca_results, scale = 0) +
  theme_bw() +
  ggtitle("ATRT Samples (z-scored high variance genes)")
ggsave(filename = file.path(plots_dir, "atrt_expression_pca.png"),
       plot = pca_plot)
