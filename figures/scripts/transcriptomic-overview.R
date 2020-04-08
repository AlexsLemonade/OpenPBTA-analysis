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

# Read in clinical data
histologies_df <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                                  guess_max = 10000)

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

gsva_scores_df <- readr::read_tsv(gsva_scores_file)
gsva_stats_df <- readr::read_tsv(gsva_stats_file)

gsva_scores_mat <- gsva_scores_df %>%
  tidyr::spread(Kids_First_Biospecimen_ID, gsea_score) %>%
  dplyr::mutate(hallmark_name = stringr::str_replace_all(
    stringr::str_remove(hallmark_name, "HALLMARK_"), "_", " ")
  ) %>%
  tibble::column_to_rownames("hallmark_name") %>%
  as.matrix()

# What gene sets are we going to include in our heatmap?
# We'll filter based on whether the ANOVA has a significant p-value
included_genesets <- gsva_stats_df %>%
  dplyr::filter(significant_anova) %>%
  dplyr::pull(hallmark_name)

# Filter to only those gene sets
gsva_scores_mat <- gsva_scores_mat[included_genesets, ]

# Let's deal with the color palette for the heatmap we will make
na_color <- divergent_palette %>%
  dplyr::filter(color_names == "na_color")

divergent_palette  <- divergent_palette %>%
  dplyr::filter(color_names != "na_color")

divergent_col_val <- seq(from = min(gsva_scores_mat),
                         to = max(gsva_scores_mat),
                         length.out = nrow(divergent_palette))

col_fun <- circlize::colorRamp2(divergent_col_val,
                                divergent_palette$hex_codes)

# Annotation bar for short histology
# First get a data frame for the annotation and ensure samples are in the
# same order
short_histology_df <- histologies_df %>%
  dplyr::filter(experimental_strategy == "RNA-Seq",
                RNA_library == "stranded") %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                short_histology) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  as.data.frame()

gsva_scores_mat <- gsva_scores_mat[, order(colnames(gsva_scores_mat))]

# The colors supplied to ComplexHeatmap::HeatmapAnnotation need to be a named
# list of named vectors
annotation_colors <- list(
  short_histology = tibble::deframe(histology_palette)
)

# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = short_histology_df,
  col = annotation_colors,
  name = "short histology",
  show_legend = FALSE
)

## Heatmap itself!
gsva_heatmap <- ComplexHeatmap::Heatmap(
  gsva_scores_mat,
  col = col_fun,
  name = "GSVA Scores",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 9),
  top_annotation = column_heatmap_annotation,
)

