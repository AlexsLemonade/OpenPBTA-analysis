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

if (!("matrixStats" %in% installed.packages())) {
  install.packages("matrixStats")
}

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
    root_dir,
    "analyses",
    "molecular-subtyping-ATRT",
    "atrt-subset",
    "atrt_zscored_expression.RDS"
  ))

# Read in final output data.frame from `01-ATRT-molecular-subtyping-data-prep.Rmd`
final_df <-
  readr::read_tsv(file.path(results_dir, "ATRT_molecular_subtypes.tsv"))

#### Plot initial heatmap ----------------------------------

# Read in initial heatmap data.frame
initial_ha_atrt_df <-
  readr::read_rds(file.path(results_dir, "initial_heatmap_annotation.RDS"))

# Define target overexpressed gene vectors
tyr_genes <-
  c("TYR",
    "MITF",
    "DCT",
    "VEGFA",
    "DNAH11",
    "SPEF1",
    "POU3F4",
    "POU3F2",
    "PBX1")
shh_genes <-
  c(
    "MYCN",
    "GLI2",
    "CDK6",
    "ASCL1",
    "HES5/6",
    "DLL1/3",
    "ZBTB7A",
    "RXF3",
    "RXF2",
    "MYBL2",
    "MXI1",
    "MEIS3",
    "MEIS2",
    "MAX",
    "INSM1",
    "FOXK1"
  )
myc_genes <-
  c(
    "MYC",
    "HOTAIR",
    "HOX",
    "TCF7L2",
    "STAT1",
    "REST",
    "RARG",
    "RAD21",
    "NR4A2",
    "IRF9",
    "IRF8",
    "FOXC1",
    "CEBPB",
    "ATF4"
  )

# Filter to only the genes of interest
zscored_expression <- zscored_expression[which(
  rownames(zscored_expression) %in% c(tyr_genes, shh_genes, myc_genes)
), ]

# Create an absolute value zscored expression matrix
zscored_expression <- as.matrix(abs(zscored_expression))

# Arrange zscored expression by rows with largest range of zscore values 
zscored_expression <- zscored_expression %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_symbol") %>%
  dplyr::mutate(row_max = matrixStats::rowMaxs(zscored_expression)) %>%
  dplyr::arrange(dplyr::desc(row_max)) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  dplyr::select(-row_max) %>%
  as.matrix()

# Filter `zscored_expression` matrix to include only the top third rows with the
# largest range of zscore values
zscored_expression <- zscored_expression[1:10,]

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
  zscored_expression,
  heatmap_legend_param = list(title = "absolute zscore value"),
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

# Create an absolute value zcored expression matrix
atrt_expression_mat <- as.matrix(abs(atrt_expression_mat))

# Arrange `atrt_expression_mat` by rows with largest range of zscore values 
atrt_expression_mat <- atrt_expression_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_symbol") %>%
  dplyr::mutate(row_max = matrixStats::rowMaxs(atrt_expression_mat)) %>%
  dplyr::arrange(dplyr::desc(row_max)) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  dplyr::select(-row_max) %>%
  as.matrix()

# Filter `atrt_expression_mat` matrix to include only the top third rows with the
# largest range of zscore values
atrt_expression_mat <- atrt_expression_mat[1:10,]

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
  atrt_expression_mat,
  heatmap_legend_param = list(title = "absolute zscore value"),
  top_annotation = final_ha_atrt
)
dev.off()
