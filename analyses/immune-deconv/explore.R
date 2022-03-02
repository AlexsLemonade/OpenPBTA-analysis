palette_file <- file.path("..", "..", "figures", "palettes", "broad_histology_cancer_group_palette.tsv")
palette_df <- read_tsv(palette_file)
divergent_file <- file.path("..", "..", "figures", "palettes", "divergent_color_palette.tsv")
divergent <- read_tsv(divergent_file)
divergent_hex <- divergent$hex_codes


as_tibble(deconv_output) %>%
  inner_join(
    select(palette_df, 
           broad_histology, 
           broad_histology_display,
           broad_histology_order, 
           cancer_group, 
           cancer_group_display)
  ) %>%
  select(-display_group, -method, -molecular_subtype) -> deconv

deconv_scores <- deconv %>%
  filter(str_detect(cell_type, "score"))
         
deconv_cells <- deconv %>%
  filter(!str_detect(cell_type, "score"))


# Find order for cancer groups, so go by BH order
palette_df %>% 
  select(broad_histology, broad_histology_order, cancer_group_display) %>% 
  filter(cancer_group_display != "Other") %>%
  arrange(broad_histology_order) %>% 
  distinct() %>% 
  mutate(cancer_group_order = 1:n()) %>%
  select(contains("cancer")) -> cancer_group_order

# Prep heatmap input
deconv_cells %>%
  drop_na(cancer_group_display, fraction) %>%
  # Only keep w/ certain sample size
  filter(cancer_group_display != "Other") %>%
  group_by(cancer_group_display, cell_type) %>%
  summarize(median_frac = median(fraction)) %>% 
  ungroup() %>%
  group_by(cancer_group_display) %>%
  # Normalize fractions as Z-scores, why not!
  mutate(zscore = (median_frac - mean(median_frac))/ sd(median_frac)) %>%
  select(-median_frac) %>%
  # get cancer order
  #inner_join(cancer_group_order) %>%
  spread(key = cancer_group_display, value = zscore) %>%
  column_to_rownames("cell_type") %>%
  t()-> for_heatmap

# ensure normalization with some tolerance
proceed <- all(colSums(for_heatmap) <= 1e-6) 
if (!(proceed)) stop("Incorrect z-scoring for heatmap")


# order the cancer groups based on their overall broad_histology groupings
for_heatmap <- for_heatmap[, cancer_group_order$cancer_group_display]

# df for heatmap annotation
annotation_df <- palette_df %>%
  select(cancer_group_display,
         broad_histology_display,
         broad_histology_order) %>%
  filter(cancer_group_display != "Other") %>%
  column_to_rownames("cancer_group_display")
annotation_df <- annotation_df[cancer_group_order$cancer_group_display, ] %>%
  select(-broad_histology_order,  # Drop extraneous column
         # Rename for display purposess
         `Broad histology` = broad_histology_display)


complex_annotation <- HeatmapAnnotation(
  df = as.data.frame(annotation_df),
  name = "Broad Histology",
  annotation_name_side = "left",
  show_legend = FALSE
)
heatmap_colors <- colorRampPalette(divergent_hex)(100)


Heatmap(
  for_heatmap,
  #col = heatmap_colors,
  name = "Immune",
  cluster_columns = FALSE,
  heatmap_legend_param = list(direction = "horizontal")
)
immune_heatmap


for_heatmap %>%
  pheatmap(fontsize = 10,
           scale = "column", angle_col = 45,
           color = heatmap_colors,
           cluster_rows=F, cluster_cols = F,
           annotation_legend = T)



par(mfrow=c(1,2))


