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
  full_join( # use full, not inner
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


as_tibble(deconv_output) %>%
  select(-display_group, -method) %>%
  filter(!str_detect(cell_type, "score")) %>%
  full_join( # use full, not inner
    select(histology_label_mapping, broad_histology, broad_histology_display, cancer_group, cancer_group_display)
  ) %>%
  ungroup() -> data

data %>% 
  drop_na() %>%
  filter(cancer_group_display != "Other") %>%
  group_by(cancer_group_display, molecular_subtype, cell_type) %>%
  #summarize(cov = sd(fraction)/mean(fraction)) %>%
  ungroup() %>%
  ggplot() + 
  aes(x = fraction) + geom_histogram() + 
  facet_grid(rows = vars(cancer_group_display), 
             cols = vars(cell_type), scales = "free")
    

# all the cell types
data %>%
  select(cell_type) %>%
  distinct() %>%
  arrange(cell_type) %>%
  pull() -> all_cell_types


# Note: there are 35 samples without a molecular subtype
data %>%
  drop_na(cancer_group_display, molecular_subtype) %>%
  filter(cancer_group_display != "Other") %>%
  inner_join(
    select(histology_label_mapping, cancer_group_display, cancer_group_hex)
  ) %>%
  # use their cell types, but NOT NK since everything is at 0 except 1 outlier. Not interesting.
  #inner_join(cell_type_names, by = "cell_type") %>%
  #filter(cell_type_name != "NK") %>%
  ## OR!! let's cycle through and find good ones
  #filter(cell_type %in% all_cell_types[1:8]) %>%
  mutate(cell_type = str_wrap(cell_type, 15)) %>%
  ggplot() + 
    aes(x = molecular_subtype, y = zscore, color = cancer_group_hex) + 
    #stat_summary(size = 0.1) + 
    geom_jitter(size = 0.25, width = 0.1) + 
  geom_hline(yintercept = 0, alpha = 0.2) +
  scale_color_identity() + 
    facet_grid(rows = vars(cell_type), scales = "free") +
  ggpubr::theme_pubr() + 
  cowplot::panel_border() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust =1),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 5)
  ) -> woah

ggsave("woah-zscore.pdf", woah, height = 25, width = 15)

