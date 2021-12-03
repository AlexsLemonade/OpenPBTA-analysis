# J. Taroni for ALSF CCDL 2020
#
# Makes a pdf panels for an RNA-seq overview figure

library(tidyverse)
library(ComplexHeatmap)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

#### PDF panels ----------------------------------------------------------------
# Here we'll set up the file names used throughout
umap_pdf <- file.path(output_dir, "umap_panel.pdf")
gsva_pdf <- file.path(output_dir, "gsva_panel.pdf")
legend_pdf <- file.path(output_dir, "broad_histology_legend.pdf")

#### Color palettes ------------------------------------------------------------

palette_dir <- file.path(root_dir, "figures", "palettes")

# Import standard color palettes for project
histology_label_mapping <- readr::read_tsv(
  file.path(palette_dir, "broad_histology_cancer_group_palette.tsv"))

# Get the broad_histology and cancer_group color mappings for the samples being
# plotted here
palette_mapping_df <- histologies_df %>%
  # Drop any rows without a broad_histology value (this is in essence dropping)
  # Normal sample and we're only working with stranded RNA-seq samples here
  filter(!is.na(broad_histology),
         experimental_strategy == "RNA-Seq",
         RNA_library == "stranded") %>%
  # Identifiers
  select(Kids_First_Biospecimen_ID,
         Kids_First_Participant_ID,
         sample_id,
         broad_histology,
         cancer_group) %>%
  # Add in hex codes & display grouping
  left_join(histology_label_mapping,
            by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Biospecimen_ID,
         broad_histology_display,
         broad_histology_hex,
         broad_histology_order) %>%
  mutate(broad_histology_display = forcats::fct_reorder(broad_histology_display,
                                                        broad_histology_order))

# Palette for GSVA scores
divergent_palette <- read_tsv(file.path(palette_dir,
                                        "divergent_color_palette.tsv"))

# The NA color is consistent across palettes
na_color <- divergent_palette %>% filter(color_names == "na_color")

# The palettes we will use for the heatmap needs NA color removed and specified
# separately
divergent_palette  <- divergent_palette %>% filter(color_names != "na_color")

#### Data prep for column heatmap annotation -----------------------------------

# We will construct an annotation bar for display_group - this is used in the
# GSVA heatmap. In this section, we prep the data frame and color palette
# required for that.

# Get a distinct version of the color keys
# This has already in essence been filtered to exclude categories that are
# irrelevant to the samples being plotted
broad_histology_color_df <- palette_mapping_df %>%
  dplyr::select(broad_histology_display, broad_histology_hex) %>%
  dplyr::distinct()

# Make color key specific to these samples
annotation_colors <- broad_histology_color_df$broad_histology_hex
names(annotation_colors) <- broad_histology_color_df$broad_histology_display

# Will no longer need this
rm(broad_histology_color_df)

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
  select(Kids_First_Biospecimen_ID, X1, X2) %>%
  # Join on color palette info
  dplyr::inner_join(palette_mapping_df,
                    by = "Kids_First_Biospecimen_ID") %>%
  select(X1, X2, broad_histology_display) %>%
  plot_dimension_reduction(point_color = "broad_histology_display",
                           x_label = "UMAP1",
                           y_label = "UMAP2",
                           alpha_value = 0.5,
                           color_palette = annotation_colors) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size = 10),
        legend.position = "none")

# Save temporary PNG
ggsave(umap_pdf, plot = umap_plot, width = 4, height = 4)

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

# We're going to drop Other from this plot and make sure it is ordered by
# broad histology display order
gsva_ordered_bsids <- palette_mapping_df %>%
  filter(broad_histology_display != "Other") %>%
  arrange(broad_histology_order) %>%
  pull(Kids_First_Biospecimen_ID)

# Use the vector of biospecimen IDs to select and order samples for GSVA display
gsva_scores_mat <- gsva_scores_mat[, gsva_ordered_bsids]

# Create a data frame just for setting up the heatmap annotation
gsva_annotation_df <- palette_mapping_df %>%
  select(Kids_First_Biospecimen_ID,
         broad_histology_display,
         broad_histology_order) %>%
  filter(Kids_First_Biospecimen_ID %in% gsva_ordered_bsids) %>%
  column_to_rownames("Kids_First_Biospecimen_ID")
gsva_annotation_df <- gsva_annotation_df[gsva_ordered_bsids, ] %>%
  select(-broad_histology_order,  # Drop extraneous column
         # Rename for display purposess
         broad_histology = broad_histology_display)

# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- HeatmapAnnotation(
  df = as.data.frame(gsva_annotation_df),
  name = "Broad Histology",
  col = list("broad_histology" = annotation_colors),
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

pdf(gsva_pdf, width = 5.75, height = 5)
draw(gsva_heatmap, heatmap_legend_side = "bottom")
dev.off()

#### display_group legend ----------------------------------------------------

# And now just the legend which will serve as the legend for the entire figure
display_group_legend <- Legend(
  labels = names(annotation_colors),
  legend_gp = gpar(fill = annotation_colors)
)

pdf(legend_pdf, width = 3, height = 6)
draw(display_group_legend)
dev.off()
