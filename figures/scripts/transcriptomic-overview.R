# J. Taroni for ALSF CCDL 2020
#
# Makes a multipanel plot for RNA-seq data display items: UMAP, GSVA scores,
# and immune deconvolution scores (xCell). This only handles stranded data right
# now.

library(tidyverse)
library(ComplexHeatmap)
library(multipanelfigure)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Final figure filename
output_png <- file.path(output_dir, "transcriptomic-overview.png")

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

#### Temporary PNG names -------------------------------------------------------
# Because we are mixing Heatmaps and ggplots in this figure, one of the best
# ways to control sizes is to create PNGs and then use multipanelfigure to
# put the final figure together.
#
# Here we'll set up the file names used throughout
umap_png <- file.path(output_dir, "temp_umap.png")
gsva_png <- file.path(output_dir, "temp_gsva.png")
deconv_png <- file.path(output_dir, "temp_deconv.png")
legend_png <- file.path(output_dir, "temp_legend.png")

#### Color palettes ------------------------------------------------------------

palette_dir <- file.path(root_dir, "figures", "palettes")

# Import standard color palettes for project
histology_label_mapping <- readr::read_tsv(
  file.path(palette_dir, "histology_label_color_table.tsv")) %>%
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>%
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

divergent_palette <- read_tsv(file.path(palette_dir,
                                        "divergent_color_palette.tsv"))
gradient_palette <- read_tsv(file.path(palette_dir,
                                       "gradient_color_palette.tsv"))

# The NA color is consistent across palettes
na_color <- divergent_palette %>% filter(color_names == "na_color")

# The palettes we will use for heatmaps need NA color removed and specified
# separately
divergent_palette  <- divergent_palette %>% filter(color_names != "na_color")
gradient_palette  <- gradient_palette %>% filter(color_names != "na_color")

#### Data prep for column heatmap annotation -----------------------------------

# We will construct an annotation bar for display_group - this is used in the
# GSVA and immune deconvolution heatmaps. In this section, we prep the data
# frame and color palette required for that.

# First get a data frame for the annotation and order samples by display_group
display_group_df <- histologies_df %>%
  filter(experimental_strategy == "RNA-Seq",
         RNA_library == "stranded") %>%
  # Join on color palette info
  dplyr::left_join(histology_label_mapping,
                    by = "Kids_First_Biospecimen_ID") %>%
  # Arrange by display_order
  arrange(display_order) %>%
  # Only keep the columns we need
  dplyr::select(Kids_First_Biospecimen_ID, display_group, hex_codes) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")


  # Get a distinct version of the color keys
  histologies_color_key_df <- display_group_df %>%
    dplyr::select(display_group, hex_codes) %>%
    dplyr::distinct()

  # Make color key specific to these samples
  annotation_colors <- unique(histologies_color_key_df$hex_codes)
  names(annotation_colors) <- unique(histologies_color_key_df$display_group)

  # Drop the hex_code column so its not made into its own annotation bar
  display_group_df <- display_group_df %>%
    dplyr::select(-hex_codes)

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

# Get the umap data set up with display_groups
umap_plot <- read_tsv(rsem_umap_file) %>%
  # Join on color palette info
  dplyr::inner_join(histology_label_mapping,
                    by = "Kids_First_Biospecimen_ID") %>%
  select(X1, X2, display_group) %>%
  plot_dimension_reduction(point_color = "display_group",
                           x_label = "UMAP1",
                           y_label = "UMAP2",
                           color_palette = annotation_colors) +
  theme(text = element_text(size = 10),
        legend.position = "none")

# Save temporary PNG
ggsave(umap_png, plot = umap_plot, width = 4, height = 4)

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
  "gsva_anova_stranded_display_group.tsv"
)

gsva_scores_df <- read_tsv(gsva_scores_file)
gsva_stats_df <- read_tsv(gsva_stats_file)

# What gene sets are we going to include in our heatmap?
# We'll filter based on whether the ANOVA has a significant p-value
included_genesets <- gsva_stats_df %>%
  filter(significant_anova) %>%
  pull(hallmark_name)

gsva_scores_mat <- gsva_scores_df %>%
  spread(Kids_First_Biospecimen_ID, gsea_score) %>%
  # Only include the significant gene sets from above
  filter(hallmark_name %in% included_genesets) %>%
  # Clean up gene set names for display
  mutate(hallmark_name = stringr::str_replace_all(
    stringr::str_remove(hallmark_name, "HALLMARK_"), "_", " ")
  ) %>%
  tibble::column_to_rownames("hallmark_name") %>%
  as.matrix()

# Let's deal with the color palette for the heatmap we will make - the values
# for the palette will be based on the scores themselves
divergent_col_val <- seq(from = min(gsva_scores_mat),
                         to = max(gsva_scores_mat),
                         length.out = nrow(divergent_palette))

gsva_col_fun <- circlize::colorRamp2(divergent_col_val,
                                     divergent_palette$hex_codes)

# Order by display_group - the biospecimens are the rownames and the data
# frame used for the annotation is already ordered by display_group
gsva_scores_mat <- gsva_scores_mat[, rownames(display_group_df)]

# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- HeatmapAnnotation(
  df = display_group_df,
  name = "display_group",
  col = list("display_group" = annotation_colors),
  na_col = na_color$hex_codes,
  annotation_name_side = "left",
  show_legend = FALSE
)

## Heatmap itself!
gsva_heatmap <- Heatmap(
  gsva_scores_mat,
  col = gsva_col_fun,
  name = "GSVA Scores",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 5),
  top_annotation = column_heatmap_annotation,
  heatmap_legend_param = list(direction = "horizontal")
)

png(gsva_png, width = 5.75, height = 5, units = "in", res = 1200)
draw(gsva_heatmap, heatmap_legend_side = "bottom")
dev.off()

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

# Palette will be specific to the values
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

# Order both xCell matrices by display_group
scores_deconv_mat <- scores_deconv_mat[, rownames(display_group_df)]
cell_type_deconv_mat <- cell_type_deconv_mat[, rownames(display_group_df)]

# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- HeatmapAnnotation(
  df = display_group_df,
  name = "display_group",
  col = list("display_group" = annotation_colors),
  na_col = na_color$hex_codes,
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Scores heatmap - use the same annotation as above
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
  show_row_dend = FALSE,
  heatmap_legend_param = list(direction = "horizontal")
)

# We can use the `%v%` operator from ComplexHeatmap to make the immune
# deconvolution panel
deconv_panel <- scores_deconv_heatmap %v% cell_type_deconv_heatmap
png(deconv_png, width = 5, height = 5, units = "in", res = 1200)
draw(deconv_panel, heatmap_legend_side = "bottom")
dev.off()

#### display_group legend ----------------------------------------------------

# And now just the legend which will serve as the legend for the entire figure
display_group_legend <- Legend(
  labels = names(annotation_colors),
  legend_gp = gpar(fill = annotation_colors)
)

png(legend_png, width = 4, height = 6, unit = "in", res = 600)
draw(display_group_legend)
dev.off()

#### Assemble multipanel figure ------------------------------------------------

transcriptomic_figure <- multi_panel_figure(columns = 7,
                                            rows = 1,
                                            width = 1200,
                                            height = 300,
                                            panel_label_type = "none")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    umap_png,
                                    col = 1:2,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    gsva_png,
                                    col = 3:4,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    deconv_png,
                                    col = 5:6,
                                    scaling = "fit")

transcriptomic_figure <- fill_panel(transcriptomic_figure,
                                    legend_png,
                                    label = NULL,
                                    scaling = "fit")

save_multi_panel_figure(transcriptomic_figure, output_png)

#### Remove temporary PNGs -----------------------------------------------------

# Now that we've constructed and saved the final output we can remove the PNG
# files for the individual panels and the legend
file.remove(c(umap_png,
              gsva_png,
              deconv_png,
              legend_png))
