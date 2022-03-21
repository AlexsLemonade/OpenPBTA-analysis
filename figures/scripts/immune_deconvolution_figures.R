# This code will be integrated into fig4-panels-gsva-umap.R when ready to go.

# source the fig4 script 
source("fig4-panels-gsva-umap.R")


# A new color mappings needed for stranded AND polya
palette_mapping_df <- histologies_df %>%
  # RNA-Seq samples only, *****BOTH****** polya and stranded
  filter(experimental_strategy == "RNA-Seq") %>%
  # Identifiers
  select(Kids_First_Biospecimen_ID,
         Kids_First_Participant_ID,
         sample_id,
         broad_histology,
         cancer_group, 
         molecular_subtype) %>%
  # Add in hex codes & display grouping
  left_join(histology_label_mapping,
            by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Biospecimen_ID,
         contains("broad_histology"), 
         contains("cancer_group"),
         molecular_subtype) 


# Analysis paths and data                                                                                                         broad_histology_order))
deconv_dir <- file.path(
  analyses_dir,
  "immune-deconv",
  "results"
)
quantiseq_file <- file.path(deconv_dir, "quantiseq_deconv-output.rds")

# Load the data
quantiseq <- read_rds(quantiseq_file)

# Define output files 
subtype_facet_file <- file.path(output_dir, "cell_types-molecular_subtypes-v1.pdf")
celltype_facet_file <- file.path(output_dir, "cell_types-molecular_subtypes-v2.pdf")
barplot_file <- file.path(output_dir, "cell_types_barplot.pdf")




# Two versions of jitter/box of fractions across cell types and molecular subtypes -----------


cancer_groups_of_interest <- c("High-grade glioma astrocytoma", "Ependymoma", "Medulloblastoma")

# find the molecular subtypes in the of interest cancer groups AND excluding unclassified, with >=3 samples
subtypes_of_interest <- palette_mapping_df %>%
  filter(cancer_group_display %in% cancer_groups_of_interest, 
         !(str_detect(molecular_subtype, "To be classified"))) %>%
  count(molecular_subtype) %>%
  filter(n >= 3) %>%
  pull(molecular_subtype)


# Filter to relevant samples 
plot_data <- quantiseq %>%
  left_join(palette_mapping_df, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  filter(molecular_subtype %in% subtypes_of_interest)


# Set shared theme
theme_set(
  ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = rel(0.6)),
          strip.text = element_text(size = rel(0.7)))
)

# Plot version 1: facet by subtype. Does not show uncharacterized cells because messes with axes.
subtype_facet <- plot_data %>%
  filter(cell_type != "uncharacterized cell") %>%
  ggplot() + 
  aes(x = cell_type, y = score, color = cancer_group_hex) + 
  geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.1) + 
  geom_jitter(width = 0.15, size = 0.75, alpha = 0.6) + 
  facet_wrap(~molecular_subtype, nrow = 3, scales = "free_y") +
  scale_color_identity() +
  labs(
    x = "Cell type", 
    y = "Estimated fraction in sample"
  )
  
# Plot version 2: facet by cell type. DOES show uncharacterized cells (but can remove this here too!)
celltype_facet <- plot_data %>%
  ggplot() + 
    aes(x = molecular_subtype, y = score, color = cancer_group_hex) + 
    geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.25) + 
    geom_jitter(width = 0.15, size = 0.75, alpha = 0.6) + 
    facet_wrap(~cell_type, nrow = 3, scales = "free_y") +
    scale_color_identity() +
    labs(
      x = "Cell type", 
      y = "Estimated fraction in sample"
    ) 


ggsave(subtype_facet_file, subtype_facet, width = 10, height = 5)
ggsave(celltype_facet_file, celltype_facet, width = 10, height = 5)
  


# Barplots of the samples, ordered in order of uncharacterized
sample_order <- plot_data %>%
  filter(cell_type == "uncharacterized cell") %>%
  arrange(score) %>% 
  pull(sample) 


barplot <- plot_data %>%
  # arrange samples
  mutate(sample = fct_relevel(sample, sample_order)) %>%
  filter(cell_type != "uncharacterized cell") %>%
  ggplot() + 
  aes(x = sample, y = score, fill = cell_type) + 
  geom_col(color = "black", size = 0.25) + 
  facet_wrap(~molecular_subtype, nrow = 3, scale = "free") +
  labs(x = "Sample", 
       y = "Cell type fraction", 
       fill = "Cell type") +
  ggsci::scale_fill_simpsons() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = rel(0.8))) 
  















