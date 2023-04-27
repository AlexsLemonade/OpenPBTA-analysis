# J. Taroni for ALSF CCDL 2020
#
# Makes a pdf panels for an RNA-seq overview figure

library(tidyverse)
library(ComplexHeatmap)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig5", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Zenodo CSV output directory and file paths
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig5a_csv <- file.path(zenodo_tables_dir, "figure-5a-data.csv")
fig5b_csv <- file.path(zenodo_tables_dir, "figure-5b-data.csv")


# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

#### PDF panels ----------------------------------------------------------------
# Here we'll set up the file names used throughout
umap_pdf <- file.path(output_dir, "umap_panel.pdf")
gsva_pdf <- file.path(output_dir, "gsva_panel.pdf")
umap_legend_pdf <- file.path(output_dir, "broad_histology_legend.pdf") ## goes with UMAP
gsva_legend_pdf <- file.path(output_dir, "cancer_group_legend.pdf")    ## does with GSVA

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
         # Keep both broad_histology for the UMAP legend, and cancer_group for the GSVA legend
         broad_histology_display,
         broad_histology_hex,
         broad_histology_order,
         cancer_group_display,
         cancer_group_hex) %>%
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
# We do this for _both_ cancer group and broad histology colors, 
# since each is used in a different plot
broad_histology_color_df <- palette_mapping_df %>%
  dplyr::select(broad_histology_display, broad_histology_hex) %>%
  dplyr::distinct() 

# Make color key specific to these samples
annotation_colors_bh<- broad_histology_color_df$broad_histology_hex
names(annotation_colors_bh) <- broad_histology_color_df$broad_histology_display

# And now for cancer group:
cancer_group_color_df <- palette_mapping_df %>%
  dplyr::select(cancer_group_display, cancer_group_hex) %>%
  dplyr::distinct() %>%
  # remove NA and "Other" displays
  tidyr::drop_na() %>%
  filter(cancer_group_display != "Other")

# Make color key specific to these samples
annotation_colors_cg<- cancer_group_color_df$cancer_group_hex
names(annotation_colors_cg) <- cancer_group_color_df$cancer_group_display

# Will no longer need these
rm(broad_histology_color_df)
rm(cancer_group_color_df)


#### UMAP plot -----------------------------------------------------------------

rsem_umap_file <- file.path(
  analyses_dir,
  "transcriptomic-dimension-reduction",
  "results",
  "rsem_stranded_log_umap_scores_aligned.tsv"
)

# Get the umap data set up with display_groups
umap_plot_df <- read_tsv(rsem_umap_file) %>%
  select(Kids_First_Biospecimen_ID, X1, X2) %>%
  # Join on color palette info
  dplyr::inner_join(palette_mapping_df,
                    by = "Kids_First_Biospecimen_ID") 


# Plot directly (without plot_dimension_reduction) to be able to specify a 
# non-mapped custom point size
umap_plot <- ggplot(umap_plot_df) + 
  aes(x = X1, y = X2, color = broad_histology_hex) + 
  geom_point(size = 2, alpha = 0.5) + 
  scale_color_identity() + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(size = 16),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )

# Save PDF
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
  "gsva_anova_stranded_cancer_group.tsv.gz"
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


# Determine cancer group counts so we can have a variable to order rows by
# But, we should keep "Other" at the _end_ regardless of counts
palette_mapping_df_ordered <- palette_mapping_df %>%
  # Remove NA and "Other" display groups
  drop_na(cancer_group_display) %>%
  filter(cancer_group_display != "Other") %>%
  # Arrange in reverse order of counts and make that the cancer_group_order
  count(cancer_group_display) %>%
  # force "Other" to be zero so we can arrange by `n` and it's at the bottom
  arrange(-n) %>%
  # Now add the proper order.
  mutate(cancer_group_order = 1:n()) %>%
  # Not needed
  select(-n) %>%
  # Join back up
  inner_join(palette_mapping_df)


# Make sure it is ordered by cancer_group order, which is based on number of samples
# Note: NA/"Other" groups are already removed above, so do not need to drop here
gsva_ordered_bsids <- palette_mapping_df_ordered %>%
  arrange(cancer_group_order) %>%
  pull(Kids_First_Biospecimen_ID)

# Use the vector of biospecimen IDs to select and order samples for GSVA display
gsva_scores_mat <- gsva_scores_mat[, gsva_ordered_bsids]

# Create a data frame just for setting up the heatmap annotation
gsva_annotation_df <- palette_mapping_df_ordered %>%
  select(Kids_First_Biospecimen_ID,
         cancer_group_display,
         cancer_group_order) %>%
  filter(Kids_First_Biospecimen_ID %in% gsva_ordered_bsids) %>%
  column_to_rownames("Kids_First_Biospecimen_ID")
gsva_annotation_df <- gsva_annotation_df[gsva_ordered_bsids, ] %>%
  select(-contains("broad_histology"),  # Drop extraneous columns
         -cancer_group_order,
         # Rename for display purposes
         cancer_group = cancer_group_display)



# Annotation bar intended for the top of the heatmap
column_heatmap_annotation <- HeatmapAnnotation(
  df = as.data.frame(gsva_annotation_df),
  name = "Cancer Group",
  col = list("cancer_group" = annotation_colors_cg),
  na_col = na_color$hex_codes,
  show_legend = FALSE,
  show_annotation_name =TRUE
)
## Heatmap itself!
gsva_heatmap <- Heatmap(
  gsva_scores_mat,
  col = gsva_col_fun,
  name = "GSVA Scores",
  na_col = na_color$hex_codes,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 5.8),
  top_annotation = column_heatmap_annotation,
  heatmap_legend_param = list(direction = "vertical")
)

pdf(gsva_pdf, width = 7, height = 4)
draw(gsva_heatmap, heatmap_legend_side = "right")
dev.off()



#### broad_histology legend ----------------------------------------------------

# UMAP legend
# note that legend text cannot be resized in `Legend()`, only the title text
# Legend() is not great at text wrapping so we're using a hack

histology_order <- names(annotation_colors_bh)

display_group_legend_plot <- tibble::as_tibble(annotation_colors_bh, 
                  rownames = "histology") %>%
  mutate(histology = factor(histology, 
                            levels = names(annotation_colors_bh))) %>%
  ggplot() + 
  aes(x = annotation_colors_bh, y = annotation_colors_bh, color = histology) + 
  geom_point(size = 3) + 
  scale_color_manual(values = annotation_colors_bh, 
                     labels = function(x) str_wrap(x, 20)) + 
  ggpubr::theme_pubr() +
  guides(color = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5), 
)
  
display_group_legend <- cowplot::get_legend(display_group_legend_plot)

cowplot::save_plot(
  umap_legend_pdf, 
  cowplot::ggdraw(display_group_legend), 
  base_width = 4.75, base_height = 1
)




#### cancer_group legend ----------------------------------------------------

# Take the annotation_colors for this from the `palette_mapping_df_ordered` tibble since it is ordered properly there
annotation_colors_cg_legend <- unique(palette_mapping_df_ordered$cancer_group_hex)
names(annotation_colors_cg_legend) <- unique(palette_mapping_df_ordered$cancer_group_display)

# GSVA legend
display_group_legend <- Legend(
  labels = names(annotation_colors_cg_legend),
  legend_gp = gpar(fill = annotation_colors_cg_legend), 
  ncol = 3,
  # reduce space between legend columns
  gap = unit(0.2, "mm")
)

pdf(gsva_legend_pdf, width = 7, height = 1.125)
draw(display_group_legend)
dev.off()



## Export CSV files for Zenodo


# 5A: UMAP
umap_plot_df %>%
  # rename UMAP columns
  dplyr::rename(UMAP1 = X1, UMAP2 = X2) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove cancer group columns since this plot is of broad histology
  dplyr::select(-contains("cancer_group")) %>%
  # export
  readr::write_csv(fig5a_csv)


# 5B: GSVA plot

# We'll transpose the gsva matrix first, because specimens are column names in the matrix
# We want them _in a column_ to arrange on
gsva_scores_mat %>%
  t() %>%
  tibble::as_tibble(rownames = "Kids_First_Biospecimen_ID") %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig5b_csv)



