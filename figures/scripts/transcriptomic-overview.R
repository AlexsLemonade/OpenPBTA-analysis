# J. Taroni for ALSF CCDL 2020
#
# Makes a multipanel plot for RNA-seq data display items: UMAP, GSVA scores,
# and immune deconvolution scores (xCell). This only handles stranded data right
# now.

library(tidyverse)
library(ComplexHeatmap)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

#### Color palettes ------------------------------------------------------------

palette_dir <- file.path(root_dir, "figures", "palettes")
histology_palette <- read_tsv(file.path(palette_dir,
                                               "histology_color_palette.tsv"))
divergent_palette <- read_tsv(file.path(palette_dir,
                                               "divergent_color_palette.tsv"))
gradient_palette <- read_tsv(file.path(palette_dir,
                                              "gradient_color_palette.tsv"))

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

umap_plot <- read_tsv(rsem_umap_file) %>%
  select(X1, X2, short_histology) %>%
  replace_na(list(short_histology = "none")) %>%
  plot_dimension_reduction(point_color = "short_histology",
                           x_label = "UMAP1",
                           y_label = "UMAP2",
                           color_palette = histology_palette)

#### GSVA scores ---------------------------------------------------------------

gsva_dir <- file.path(
  analyses_dir,
  "gene-set-enrichment-analysis"
)

gsva_scores_file <- file.path(
  gsva_dir,
  "results",
  "gsva_scores_stranded.tsv"
)

gsva_stats_file <- file.path(
  gsva_dir,
  "results",
  "gsva_anova_stranded_short_histology.tsv"
)

gsva_scores_df <- read_tsv(gsva_scores_file)
gsva_stats_df <- read_tsv(gsva_stats_file)

gsva_scores_mat <- gsva_scores_df %>%
  spread(Kids_First_Biospecimen_ID, gsea_score) %>%
  mutate(hallmark_name = stringr::str_replace_all(
    stringr::str_remove(hallmark_name, "HALLMARK_"), "_", " ")
  ) %>%
  tibble::column_to_rownames("hallmark_name") %>%
  as.matrix()

# What gene sets are we going to include in our heatmap?
# We'll filter based on whether the ANOVA has a significant p-value
included_genesets <- gsva_stats_df %>%
  filter(significant_anova) %>%
  pull(hallmark_name)

# Filter to only those gene sets
gsva_scores_mat <- gsva_scores_mat[included_genesets, ]

# Let's deal with the color palette for the heatmap we will make
na_color <- divergent_palette %>%
  filter(color_names == "na_color")

divergent_palette  <- divergent_palette %>%
  filter(color_names != "na_color")

divergent_col_val <- seq(from = min(gsva_scores_mat),
                         to = max(gsva_scores_mat),
                         length.out = nrow(divergent_palette))

gsva_col_fun <- circlize::colorRamp2(divergent_col_val,
                                     divergent_palette$hex_codes)

# Annotation bar for short histology
# First get a data frame for the annotation and ensure samples are in the
# same order
short_histology_df <- histologies_df %>%
  filter(experimental_strategy == "RNA-Seq",
                RNA_library == "stranded") %>%
  select(Kids_First_Biospecimen_ID,
                short_histology) %>%
  arrange(Kids_First_Biospecimen_ID) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  as.data.frame()

gsva_scores_mat <- gsva_scores_mat[, order(colnames(gsva_scores_mat))]

# The colors supplied to HeatmapAnnotation need to be a named
# list of named vectors
annotation_colors <- list(
  short_histology = tibble::deframe(histology_palette)
)

# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- HeatmapAnnotation(
  df = short_histology_df,
  col = annotation_colors,
  name = "short histology",
  show_legend = FALSE
)

## Heatmap itself!
gsva_heatmap <- Heatmap(
  gsva_scores_mat,
  col = gsva_col_fun,
  name = "GSVA Scores",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 9),
  top_annotation = column_heatmap_annotation
)

#### Immune deconvolution ------------------------------------------------------

immune_deconv_dir <- file.path(
  analyses_dir,
  "immune-deconv"
)

immune_deconv_file <- file.path(
  immune_deconv_dir,
  "results",
  "deconv-output-for-figures.RData"
)

load(immune_deconv_file)

# Get in wide format for making a heatmap
deconv_mat <- deconv.res %>%
  select(cell_type, `sample`, fraction) %>%
  spread(`sample`, fraction) %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

# We only want the stranded data - this particular module combined the two
# RNA-seq datasets
deconv_mat <- deconv_mat[, colnames(gsva_scores_mat)]

# We need a new color palette, etc.
na_color <- gradient_palette %>%
  filter(color_names == "na_color")

gradient_palette  <- gradient_palette %>%
  filter(color_names != "na_color")

gradient_col_val <- seq(from = min(deconv_mat),
                         to = max(deconv_mat),
                         length.out = nrow(gradient_palette))

deconv_col_fun <- circlize::colorRamp2(gradient_col_val,
                                       gradient_palette$hex_codes)

# There are 3 general scores from xCell - immune score, microenvironment score,
# and stroma score
# We're going to separate those rows out from the cell type scores
scores_deconv_mat <- deconv_mat[grep("score", rownames(deconv_mat)), ]
cell_type_deconv_mat <- deconv_mat[-grep("score", rownames(deconv_mat)), ]

# To keep ordering consistent across the two xCell heatmaps, we'll order by
# histology
short_histology_df <- short_histology_df %>%
  # We lose the rownames if we use arrange without this step
  tibble::rownames_to_column("biospecimen") %>%
  arrange(short_histology)

# Order both xCell matrices by histology
scores_deconv_mat <- scores_deconv_mat[, short_histology_df$biospecimen]
cell_type_deconv_mat <- cell_type_deconv_mat[, short_histology_df$biospecimen]

# Make the biospecimen ID the rownames again
short_histology_df <- short_histology_df %>%
  tibble::column_to_rownames("biospecimen")

# Annotation bar intended for the top of the heatmap needs to be remade because
# of the ordering but the colors are the same
column_heatmap_annotation <- HeatmapAnnotation(
  df = short_histology_df,
  col = annotation_colors,
  name = "short histology",
  show_legend = FALSE
)

# Scores heatmap
scores_deconv_heatmap <- Heatmap(
  scores_deconv_mat,
  col = deconv_col_fun,
  name = "xCell fraction",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 9),
  top_annotation = column_heatmap_annotation,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = FALSE
)

# For display filter to cell types with highest variance
row_variances <- matrixStats::rowVars(cell_type_deconv_mat)
high_var_index <- which(row_variances > quantile(row_variances, 0.75))
cell_type_deconv_mat <- cell_type_deconv_mat[high_var_index, ]

# Cell type heatmap
cell_type_deconv_heatmap <- Heatmap(
  cell_type_deconv_mat,
  col = deconv_col_fun,
  name = "xCell fraction",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 9),
  cluster_columns = FALSE,
  show_row_dend = FALSE
)

# We can use the `%v%` operator from ComplexHeatmap to make the immune
# deconvolution panel
deconv_panel <- scores_deconv_heatmap %v% cell_type_deconv_heatmap

#### PNGs ----------------------------------------------------------------------
# We're mixing Heatmaps and ggplots - we can write everything to PNG and then
# reassemble because viewports were not playing nicely :(







