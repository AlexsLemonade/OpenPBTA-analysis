# Oncoprint Landscape Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 3.

#### Set Up --------------------------------------------------------------------

# Install maftools
if (!("maftools" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("maftools")
}
library(maftools)

# Load in patchwork for assembling the final multipanel figure
library(patchwork)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to the data obtained via `00-map-to-sample_id.R`
scratch_oncoprint_dir <- file.path(root_dir, "scratch", "oncoprint_files")

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Path to oncoprint landscape directory
oncoprint_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Source the color palette for the oncoprint plots
source(
  file.path(
    oncoprint_dir,
    "util",
    "oncoplot-palette.R"
  )
)

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF files
maf_df_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_maf.tsv"
  ))

maf_df_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_maf.tsv"
  ))

# Read in cnv files
cnv_file_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_cnv.tsv"
  ))

cnv_file_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_cnv.tsv"
  ))

# Read in fusion files
fusion_file_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_fusions.tsv"
  ))

fusion_file_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_fusions.tsv"
  ))

# Read in gene list
gene_list <-
  readr::read_tsv(
    file.path(
      root_dir,
      "analyses",
      "interaction-plots",
      "results",
      "gene_disease_top50.tsv"
    )
  ) %>%
  dplyr::pull("gene")

#### Functions ----------------------------------------------------------

source(file.path(
  oncoprint_dir,
  "util",
  "oncoplot-functions.R"
))

#### Generate Oncoprints ------------------------------------------------------

# Generate the primary-plus oncoprint
png(
  file.path(output_dir, "fig3a-oncoprint-lanscape.png"),
  width = 65,
  height = 30,
  units = "cm",
  res = 300
)
prepare_and_plot_oncoprint(
  maf_df_primary_plus,
  cnv_file_primary_plus,
  metadata,
  fusion_file_primary_plus,
  goi_list = gene_list,
  color_palette = color_palette
)
dev.off()

# Generate the primary-only oncoprint
png(
  file.path(output_dir, "fig3b-oncoprint-lanscape.png"),
  width = 65,
  height = 30,
  units = "cm",
  res = 300
)
prepare_and_plot_oncoprint(
  maf_df_primary_only,
  cnv_file_primary_only,
  metadata,
  fusion_file_primary_only,
  goi_list = gene_list,
  color_palette = color_palette
)
dev.off()

#### Assemble multipanel plot -------------------------------------------------

# TODO: Determine the best way to do this as pathwork::plot_layout was not
#       successful
