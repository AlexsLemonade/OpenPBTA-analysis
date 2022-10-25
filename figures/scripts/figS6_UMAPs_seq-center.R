# S. Spielman 2022
# Script to make UMAPs exploring potential sequencing center batch effects

# Libraries --------
library(magrittr)

# Define paths and files ---------------------
dim_red_dir <- file.path(
  "analyses",
  "transcriptomic-dimension-reduction"
)

metadata_file <- file.path("data", "pbta-histologies.tsv")
pal_file <- file.path("figures", "palettes", "broad_histology_cancer_group_palette.tsv")
rsem_umap_file <- file.path(
  dim_red_dir,
  "results",
  "rsem_stranded_log_umap_scores_aligned.tsv"
)


# Output:
output_dir <- file.path("figures", "pdfs", "supp", "figS5", "panels")
umap_all_bh_pdf <- file.path(output_dir, "umap_sequencing_all_histologies.pdf")
umap_shared_bh_pdf <- file.path(output_dir, "umap_sequencing_shared_histologies.pdf")
umap_centers_only_pdf <- file.path(output_dir, "umap_sequencing_centers.pdf")


# Source the UMAP plotting function -------------
source(file.path(dim_red_dir, "util", "dimension-reduction-functions.R"))


# Read in data files ----------------
histologies_df <- readr::read_tsv(metadata_file, guess_max = 100000)

rsem_umaps <- readr::read_tsv(rsem_umap_file)

palette_mapping_df <- readr::read_tsv(pal_file) %>%
  dplyr::select(broad_histology, broad_histology_display, broad_histology_hex) 


# Define the color key -----------------
broad_histology_color_df <- palette_mapping_df %>%
  dplyr::select(broad_histology_display, broad_histology_hex) %>%
  dplyr::distinct() 
annotation_colors_bh<- broad_histology_color_df$broad_histology_hex
names(annotation_colors_bh) <- broad_histology_color_df$broad_histology_display


# Render UMAPs ----------------------

# First, link IDs to histology display
palette_mapping_df <- palette_mapping_df %>%
  dplyr::inner_join(histologies_df, by="broad_histology") %>%
  dplyr::select(Kids_First_Biospecimen_ID, broad_histology_display)

# Prepare data for UMAP plot
umap_df <- rsem_umaps %>%
  dplyr::select(Kids_First_Biospecimen_ID, seq_center, X1, X2) %>%
  # Join on color palette info
  dplyr::left_join(palette_mapping_df,
                    by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(X1, X2, seq_center, broad_histology_display) 

# Shuffle the rows going in for alpha purposes, and
#  set order for `seq_center` as group size
set.seed(2022)
umap_df <- umap_df %>%
  dplyr::sample_n(nrow(umap_df)) %>%
  dplyr::mutate(seq_center = forcats::fct_infreq(seq_center))

# UMAP of all histologies
umap_all_bh <- plot_dimension_reduction(umap_df,
                                        point_color = "broad_histology_display",
                                        point_shape = "seq_center",
                                        x_label = "UMAP1",
                                        y_label = "UMAP2",
                                        alpha_value = 0.15) + 
  ggplot2::labs(title = "All broad histologies", 
                color = "Broad histology", 
                shape = "Sequencing center")

# UMAP of all _shared_ histologies

# first, which histologies are shared, aka processed at >1 center?
shared_histologies <- umap_df %>% 
  dplyr::count(seq_center, broad_histology_display) %>%
  dplyr::count(broad_histology_display) %>%
  dplyr::filter(n>1) %>%
  dplyr::pull(broad_histology_display)

umap_shared_bh <- plot_dimension_reduction(dplyr::filter(umap_df, broad_histology_display %in% shared_histologies),
                                           point_color = "broad_histology_display",
                                           point_shape = "seq_center",
                                           x_label = "UMAP1",
                                           y_label = "UMAP2",
                                           alpha_value = 0.3) +
  ggplot2::labs(title = "Broad histologies sequenced at multiple centers", 
                color = "Broad histology", 
                shape = "Sequencing center")

# UMAP colored by sequencing center, ignoring histology:
umap_centers_only <- plot_dimension_reduction(umap_df,
                                              point_color = "seq_center",
                                              x_label = "UMAP1",
                                              y_label = "UMAP2",
                                              alpha_value = 0.2) +
  ggplot2::labs(title = "Broad histologies sequenced at multiple centers", 
                color = "Broad histology", 
                shape = "Sequencing center") +
  # colorblind friendly palette:
  colorblindr::scale_color_OkabeIto()


# Export panels ----------
ggplot2::ggsave(umap_all_bh_pdf, umap_all_bh, width = 8, height = 5)
ggplot2::ggsave(umap_shared_bh_pdf, umap_shared_bh, width = 8, height = 5)
ggplot2::ggsave(umap_centers_only_pdf, umap_centers_only, width = 8, height = 5)

