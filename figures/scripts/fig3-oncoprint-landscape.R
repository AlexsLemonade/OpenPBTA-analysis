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

# Path to the data obtained via `bash download-data.sh`.
data_dir <- file.path(root_dir, "data")

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Path to oncoprint landscape directory
oncoprint_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Source the color palette for plots
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
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_maf.tsv"
    )
  )

maf_df_primary_only <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_maf.tsv"
    )
  )

# Read in cnv files
cnv_file_primary_plus <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_cnv.tsv"
    )
  )

cnv_file_primary_only <-
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_cnv.tsv"
    )
  )

# Read in fusion files
fusion_file_primary_plus <- 
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-plus_fusions.tsv"
    )
  )

fusion_file_primary_only <- 
  readr::read_tsv(
    file.path(
      root_dir,
      "scratch",
      "oncoprint_files",
      "all_participants_primary-only_fusions.tsv"
    )
  )

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

#### Format the color palette -----------------------------------------------------------

# Read in the histology color palette from the `figures/palettes` directory
histology_col_palette <- readr::read_tsv(
  file.path(root_dir, "figures", "palettes", "histology_color_palette.tsv")) %>%
  as.data.frame() %>%
  # Store the histology names as row names for use with the `oncoplot`` function
  tibble::column_to_rownames("color_names")

# Make the sourced color palette a data.frame and join the `histology_col_palette`
# values -- this paletter will have hex codes for short histologies, SNVs, CNVs,
# and fusion data categories
color_palette <- as.data.frame(color_palette) %>%
  select(hex_codes = color_palette) %>%
  rbind(col_palette)

#### Generate Oncoprints ------------------------------------------------------

primary_plus_oncoprint <- prepare_and_plot_oncoprint(maf_df_primary_plus,
                                                     cnv_file_primary_plus,
                                                     fusion_file_primary_plus,
                                                     gene_list = gene_list,
                                                     color_palette = color_palette)

primary_only_oncoprint <- prepare_and_plot_oncoprint(maf_df_primary_only,
                                                     cnv_file_primary_only,
                                                     fusion_file_primary_only,
                                                     gene_list = gene_list,
                                                     color_palette = color_palette)

#### Assemble multipanel plot -------------------------------------------------

# Combine plots with patchwork
# Layout of the two plots will be one over the other (1 column and 2 rows)
combined_plot <- primary_plus_oncoprint + primary_only_oncoprint +
  plot_layout(ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(# add uniform labels
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9))

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig3-oncoprint-lanscape.png"),
                width = 12, height = 8,
                units = "in"
)
