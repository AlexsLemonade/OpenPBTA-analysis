library(tidyverse)

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



palette_dir <- file.path(root_dir, "figures", "palettes")

# Import standard color palettes for project
histology_label_mapping <- readr::read_tsv(
  file.path(palette_dir, "broad_histology_cancer_group_palette.tsv"))


# A new color mappings for stranded and polya
palette_mapping_df <- histologies_df %>%
  # RNA-Seq samples only, *****BOTH****** polya and stranded
  filter(experimental_strategy == "RNA-Seq") %>%
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

# Analysis paths and data                                                                                                         broad_histology_order))
deconv_dir <- file.path(
  analyses_dir,
  "immune-deconv",
  "results"
)

# Data - need to load, and creates variable `deconv_output`
deconv_output_file <- file.path(deconv_dir, "deconv-output.RData")
load(deconv_output_file)



# This is a plot of the microenvironment scores, which are roughly expected to correlate NEGATIVELY with purity
as_tibble(deconv_output) %>%
  select(-display_group, -method) %>%
  filter(str_detect(cell_type, "score"), 
        (str_detect(cell_type, "micro"))) %>%
  #mutate(score = str_replace(cell_type, " score", "")) %>%
  select(-cell_type) %>%
  inner_join( 
    select(histology_label_mapping, broad_histology, broad_histology_display, cancer_group, cancer_group_display)
  ) %>%
  drop_na() %>%
  filter(cancer_group_display != "Other") %>%
  ggplot(aes(x = fct_reorder(cancer_group_display, fraction, .desc=T), 
             y = fraction)) + 
  geom_violin() + 
  stat_summary() 

# here are the microenvironment scores also faceted by subtype, jitter
as_tibble(deconv_output) %>%
  select(-display_group, -method) %>%
  filter(str_detect(cell_type, "score"), 
         (str_detect(cell_type, "micro"))) %>%
  #mutate(score = str_replace(cell_type, " score", "")) %>%
  select(-cell_type) %>%
  full_join( # use full, not inner
    select(histology_label_mapping, broad_histology, broad_histology_display, cancer_group, cancer_group_display)
  ) %>%
  drop_na() %>%
  filter(cancer_group_display != "Other") %>%
  ggplot(aes(x = fct_reorder(cancer_group_display, fraction, .desc=T), 
             y = fraction)) + 
  geom_jitter() + 
  facet_wrap(~molecular_subtype, scales = "free")



# These are the cell types used in the benchmark paper - 
# https://github.com/icbi-lab/immune_deconvolution_benchmark
# https://github.com/icbi-lab/immune_deconvolution_benchmark/blob/master/notebooks/70_pub_figures.Rmd
cell_type_names = tibble(
  cell_type      = c("B cell", "Dendritic cell", "Macrophage/Monocyte",
                     "NK cell", "T cell CD4+", "T cell CD8+", "T cell CD4+ (non-regulatory)",  "T cell regulatory (Tregs)",
                     "Monocyte", "T cell", "Cancer associated fibroblast", "Endothelial cell", "cancer cell"),
  cell_type_name = c("B",      "DC",             "Mac/Mono",
                     "NK",      "T CD4+",      "T CD8+" ,     "T CD4+ n.r.",                   "T reg",
                     "Mono",     "T",      "CAF",                          "Endo",             "Cancer")
)


hgat_epn_cancer_groups <- c("Ependymoma", "High-grade glioma astrocytoma", "Diffuse midline glioma", "Diffuse intrinsic pontine glioma")


as_tibble(deconv_output) %>%
  select(-display_group, -method) %>%
  #filter(!str_detect(cell_type, "score")) %>%
  filter(cell_type == "immune score" | str_detect(cell_type, "CD8+")) %>%
  rename(Kids_First_Biospecimen_ID = sample) %>%
  inner_join(histologies_df) %>%
  ungroup() %>%
  filter(cancer_group %in% hgat_epn_cancer_groups) %>%
  drop_na(molecular_subtype) %>%
  select(cell_type, molecular_subtype, Kids_First_Biospecimen_ID, cancer_group, fraction) %>%
  inner_join(
    select(histology_label_mapping, cancer_group, cancer_group_display, cancer_group_hex)
  ) %>%
  # use their cell types, but NOT NK since everything is at 0 except 1 outlier. Not interesting.
  #inner_join(cell_type_names, by = "cell_type") %>%
  #filter(cell_type_name != "NK") %>%
  mutate(cell_type = str_wrap(cell_type, 15)) %>%
  #filter(str_detect(cell_type, "T cell CD8+")) %>%
  ggplot() + 
    aes(x = molecular_subtype, y = fraction, color = cancer_group_hex) + 
    #stat_summary() + 
    geom_boxplot(outlier.shape=NA, color = "black", size = 0.25) + 
    geom_jitter(size = rel(0.75), width = 0.1) + 
  scale_color_identity() + 
  facet_grid(rows = vars(cell_type),
             #cols = vars(molecular_subtype),
             scales = "free") +
  ggpubr::theme_pubr() + 
  cowplot::panel_border() +
  theme(
    axis.text.x = element_text(size = 6, angle= 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7)
  ) -> plot_option


ggsave(file.path(output_dir, "temporary_box-jitter.pdf"), 
       plot_option,
       width = 10, height = 6)

