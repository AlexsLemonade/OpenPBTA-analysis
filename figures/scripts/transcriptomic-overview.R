# J. Taroni for ALSF CCDL 2020
#
# Makes a multipanel plot for RNA-seq data display items: UMAP, GSVA scores,
# and immune deconvolution scores (xCell). This only handles stranded data right
# now.

`%>%` <- dplyr::`%>%`

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

#### Color palettes ------------------------------------------------------------

palette_dir <- file.path(root_dir, "figures", "palettes")
histology_palette <- readr::read_tsv(file.path(palette_dir,
                                               "histology_color_palette.tsv"))
divergent_palette <- readr::read_tsv(file.path(palette_dir,
                                               "divergent_color_palette.tsv"))

#### UMAP plot -----------------------------------------------------------------

dim_red_dir <- file.path(
  analyses_dir,
  "transcriptomic-dimension-reduction"
)

source(file.path(dim_red_dir, "util", "dimension-reduction-functions.R"))

rsem_umap_file <- file.path(
  dim_red_dir,
  "results",
  "rsem_stranded_log_umap_scores_aligned.tsv"
)

umap_plot <- readr::read_tsv(rsem_umap_file) %>%
  dplyr::select(X1, X2, short_histology) %>%
  tidyr::replace_na(list(short_histology = "none")) %>%
  plot_dimension_reduction(point_color = "short_histology",
                           x_label = "UMAP1",
                           y_label = "UMAP2",
                           color_palette = histology_palette)

#### GSVA scores ---------------------------------------------------------------
