# S. Spielman for ALSF CCDL, 2022
#
# Makes PDF panels from the `immune-deconv` module for Figures 5 and S6, specifically:
## Panel 5C, immune cell fractions across immune cells, faceted by cancer group
## Panel 5E, CD274 expression across molecular subtypes
## Panel S6E, immune cell fractions across molecular subtypes, faceted by immune cells
## Panel S6F, CD8/CD4 ratio across a subset of molecular subtypes
# 
# All code is adapted from `analyses/immune-deconv/02-visualize_quantiseq.Rmd`


library(tidyverse)

# set seed for jitter plot reproducibility
set.seed(2022)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directories
output_dir_fig5 <- file.path(root_dir, "figures", "pdfs", "fig5", "panels")
if (!dir.exists(output_dir_fig5)) {
  dir.create(output_dir_fig5, recursive = TRUE)
}
output_dir_figS6 <- file.path(root_dir, "figures", "pdfs", "supp", "figs6", "panels")
if (!dir.exists(output_dir_figS6)) {
  dir.create(output_dir_figS6, recursive = TRUE)
}

# Input directories
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "immune-deconv", "results")
palette_dir <- file.path(root_dir, "figures", "palettes")

# Array of groups we want to visualize in 5C
cancer_groups_of_interest <- c("Pilocytic astrocytoma", 
                               "Diffuse midline glioma",
                               "Craniopharyngioma", 
                               "Ganglioglioma",
                               "Ependymoma", 
                               "Medulloblastoma",
                               "Schwannoma",
                               "Neurofibroma Plexiform",
                               "Other embryonal tumor",
                               "Atypical Teratoid Rhabdoid Tumor",
                               "Meningioma",
                               "Dysembryoplastic neuroepithelial tumor")

## Input and output files -------------------------------------

hist_file <- file.path(data_dir, "pbta-histologies.tsv")
palette_file <- file.path(palette_dir, "broad_histology_cancer_group_palette.tsv")
quantiseq_file <- file.path(analysis_dir, "quantiseq_deconv-output.rds")
polya_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")
stranded_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

# Figure 5 output:
cancer_group_pdf <- file.path(output_dir_fig5, "quantiseq-cell_types-cancer_groups.pdf")
cd274_mb_plot_pdf <- file.path(output_dir_fig5, "cd274_expression_mb_subtypes.pdf")

# Figure S6 output:
subtypes_pdf <- file.path(output_dir_figS6, "quantiseq-cell_types-molecular_subtypes.pdf")
cd8_cd4_ratio_pdf <- file.path(output_dir_figS6, "cd8_cd4_ratio.pdf")


# Prepare data for plots --------------

# Read in clinical data and palette file
histologies_df <- read_tsv(hist_file, guess_max = 10000)
palette_df <- readr::read_tsv(palette_file)

# Read in analysis data
quantiseq <- read_rds(quantiseq_file)
polya <- read_rds(polya_file)
stranded <- read_rds(stranded_file)

# Prepare palette/histology data data
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
  left_join(palette_df,
            by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Biospecimen_ID,
         contains("broad_histology"), 
         contains("cancer_group"),
         molecular_subtype)


#### Prepare quanTIseq data --> `quantiseq_subset`
# spread and create cd8/cd4 ratio
quantiseq_spread <- quantiseq %>%
  # get a column for each cell type
  spread(cell_type, score) %>%
  mutate(cd8_cd4_ratio = `T cell CD8+` / `T cell CD4+ (non-regulatory)`)

quantiseq_gather <- quantiseq_spread %>%  
  gather(cell_type, score, -c(sample, library, method))

# First, find the molecular subtypes in the of interest cancer groups AND excluding unclassified, with >=3 samples
subtypes_of_interest <- palette_mapping_df %>%
  filter(cancer_group_display %in% cancer_groups_of_interest, 
         !(str_detect(molecular_subtype, "To be classified"))) %>%
  count(molecular_subtype) %>%
  filter(n >= 3) %>%
  pull(molecular_subtype)


# Now, filter to relevant samples and remove uncharacterized fractions
quantiseq_subset <- quantiseq_gather %>%
  left_join(palette_mapping_df, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  filter(molecular_subtype %in% subtypes_of_interest, 
         cell_type != "uncharacterized cell") %>%
  # Change loss --> loss here so inherited by all
  mutate(molecular_subtype = ifelse(molecular_subtype == "DMG, H3 K28, TP53 loss", 
                                    "DMG, H3 K28, TP53 lost", 
                                    molecular_subtype))


#### Prepare expression data to explore PDL1 aka CD274 ---> `expression_pdl1`

# Combine polya and stranded to get expression for pdl1, while keeping a library annotation
polya_expression <- polya %>%
  rownames_to_column("gene") %>%
  filter(gene == "CD274") %>%
  gather(-gene, key = sample, value = expression) %>%
  mutate(library = "polya")

expression_pdl1 <- stranded %>%
  rownames_to_column("gene") %>%
  filter(gene == "CD274") %>%
  gather(-gene, key = sample, value = expression) %>%
  mutate(library = "stranded") %>%
  bind_rows(polya_expression) %>%
  as_tibble() %>%
  mutate(log2_expression = log(expression +1, 2)) %>%
  select(-expression)


## Plot 5C ----------------------

# filter to relevant samples and remove uncharacterized fractions
quantiseq_subset_cg <- quantiseq_gather %>%
  filter(cell_type != "uncharacterized cell" & cell_type != "cd8_cd4_ratio") %>%
  left_join(palette_mapping_df, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  filter(cancer_group_display %in% cancer_groups_of_interest)

# We need to create a new scheme for labeling that shows wrapped cancer groups with `(n=X)`
quantiseq_subset_cg <- quantiseq_subset_cg %>%
  count(cancer_group_display, cell_type) %>%
  select(-cell_type) %>%
  distinct() %>%
  inner_join(
    select(palette_df, contains("cancer_group"))
  ) %>%
  select(cancer_group_display, n, cancer_group_hex) %>%
  # Create wrapped with (n=X) factor column for cancer groups
  mutate(cancer_group_display_n = stringr::str_wrap(glue::glue("{cancer_group_display} (N={n})"), 30),
                cancer_group_display_n = forcats::fct_reorder(cancer_group_display_n, n, .desc=T)) %>%
  inner_join(quantiseq_subset_cg)

# plot by cancer group
cancer_group_plot <- quantiseq_subset_cg %>%
  ggplot() + 
  aes(x = cell_type, y = score, color = cancer_group_hex) + 
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, 
              # This shape is helpful when compiling
              shape = 16) + 
  geom_boxplot(outlier.size = 0, 
               size = 0.2,
               color = "black", 
               alpha = 0,
               # remove whiskers
               coef = 0) +
  facet_wrap(~cancer_group_display_n, nrow = 2, scales = "free_y") +
  scale_color_identity() + 
  labs(
    x = "Immune cell", 
    y = "Estimated fraction in sample"
  ) +
  ggpubr::theme_pubr() + 
  cowplot::panel_border() +
  theme(
    # Sizing for compilation
    axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5.5),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 4.3),
    strip.background = element_rect(size = 0.4),
    axis.line = element_line(size = rel(0.4)),
    axis.ticks = element_line(size = rel(0.4)),
    legend.position = "none"
  )

# Sized for compilation
ggsave(cancer_group_pdf, cancer_group_plot, 
       width = 8, height = 2.65, useDingbats=FALSE)


## Plot 5E ------------------------------------


# Prepare this plot data
pdl1_cd8_mb <- expression_pdl1 %>%
  inner_join(quantiseq_subset) %>%
  select(log2_expression, cell_type, library, score, molecular_subtype, cancer_group_display, cancer_group_hex) %>%
  filter(cell_type == "T cell CD8+", 
         cancer_group_display != "Other", 
         cancer_group_display == "Medulloblastoma") %>%
  # Count group sizes to use in x-axis labels
  group_by(molecular_subtype) %>%
  mutate(subtype_count = n()) %>%
  ungroup() %>%
  mutate(molecular_subtype = glue::glue("{molecular_subtype}\n(N = {subtype_count})"))

wilcox_df <- ggpubr::compare_means(log2_expression ~ molecular_subtype, 
                                   data = pdl1_cd8_mb,
                                   method = "wilcox.test") %>%
  # Add y-axis positions for P-values
  mutate(y.position = c(0.8,
                        NA,
                        1.8,
                        1.0,
                        2.0,
                        2.2))

# distributions of PDL1 expression
cd274_mb_plot <- ggplot(pdl1_cd8_mb) + 
  aes(x = molecular_subtype,
      y = log2_expression) +
  # remove outliers
  geom_boxplot(outlier.shape = NA, 
               color = "grey40", 
               size = 0.2) + 
  geom_jitter(width = 0.15, 
              size = 0.55, 
              alpha = 0.5,
              # Helpful shape for compiling
              shape = 16) + 
  ggpubr::stat_pvalue_manual(wilcox_df, 
                             label = "p = {p.adj}",
                             label.size = 1.5,
                             bracket.size = 0.125) +
  scale_color_identity() +
  labs(x = "Molecular subtype of sample",
       y = "CD274 log2(FPKM+1)") +
  ggpubr::theme_pubr() + 
  # Set sizing for compilation
  theme(
    axis.text.x = element_text(hjust = 0.55,
                               size = 4),
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 5.5),
    axis.line  = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2)
    ) 

# Panel sized for compilation
ggsave(cd274_mb_plot_pdf, 
       cd274_mb_plot, width = 2, height = 2, 
       useDingbats=FALSE)



## Plot S6E -----------------------------------

# First, remove cd8_cd4_ratio
quantiseq_no_ratio <- quantiseq_subset %>%
  filter(cell_type != "cd8_cd4_ratio")

# add an ordering to molecular subtypes based on broad_histology_display labels
quantiseq_no_ratio <- quantiseq_no_ratio %>%
  select(broad_histology_display, molecular_subtype) %>%
  unique() %>%
  arrange(broad_histology_display, molecular_subtype) %>%
  mutate(mol_subtype_order = 1:n()) %>%
  #join back to df
  inner_join(quantiseq_no_ratio, by = c("broad_histology_display", "molecular_subtype")) %>%
  dplyr::mutate(molecular_subtype = forcats::fct_reorder(molecular_subtype, mol_subtype_order)) 

# Faceted plot of subtypes, cell types. Jitter points are colored by underlying cancer group.
subtypes_celltypes_plot <- ggplot(quantiseq_no_ratio) +
  aes(x = molecular_subtype, y = score) + 
  geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.2) + 
  geom_jitter(width = 0.15, size = 0.66, alpha = 0.6, aes(color = broad_histology_hex),
              # Helpful shape for compiling
              shape = 16) + 
  facet_wrap(~cell_type, ncol = 5, scales = "free_y") +
  scale_color_identity() +
  labs(
    x = "Molecular subtype of tumor sample", 
    y = "Estimated fraction in sample"
  ) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        axis.title = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(0.7)),
        strip.text = element_text(size = rel(0.7)),
        axis.line = element_line(size = rel(0.7)),
        axis.ticks = element_line(size = rel(0.7)))

# Add hack legend 
cg_color_df <- quantiseq_no_ratio %>%
  select(broad_histology_hex, broad_histology_display) %>%
  distinct()

cg_color_list <- cg_color_df$broad_histology_hex
names(cg_color_list) <- cg_color_df$broad_histology_display

legend_plot <- ggplot(quantiseq_no_ratio) + 
  aes(x = score, y = score, color = broad_histology_display) + 
  # Set alpha to match plot point alpha
  geom_point(alpha = 0.6) +
  scale_color_manual(name = "Broad Histology", values = cg_color_list) + 
  ggpubr::theme_pubr() +
  # Need to make text very small
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  )

legend_panel <- cowplot::get_legend(legend_plot)

subtypes_plot <- cowplot::plot_grid(subtypes_celltypes_plot, legend_panel, 
                                nrow = 2, 
                                rel_heights = c(1, 0.1))

ggsave(subtypes_pdf, subtypes_plot, 
       width = 11, height = 3.5, useDingbats=FALSE)




## Plot S6F -----------------------------------


# subset for just cd8/cd4 ratio
ratio_df <- quantiseq_subset %>%
  filter(cell_type == "cd8_cd4_ratio",
         cancer_group != "Other",
         !is.nan(score),
         !is.infinite(score)
  )

# distributions of PDL1 expression
cd8_cd4_ratio_plot <- ggplot(ratio_df) + 
  aes(x = molecular_subtype,
      y = score) +
  # remove outliers
  geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.2) + 
  geom_jitter(width = 0.15, size = 0.75, alpha = 0.6,  aes(color = broad_histology_hex),
              # Helpful shape for compiling
              shape = 16) + 
  scale_color_identity() +
  labs(x = "Molecular subtype of tumor sample",
       y = "Ratio of CD8+/CD4+ T cell fractions") +
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(hjust = 1, 
                                   size = rel(0.3), 
                                   angle = 45),
        axis.text.y = element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7)),
        axis.line = element_line(size = rel(0.7)),
        axis.ticks = element_line(size = rel(0.7)))


ggsave(cd8_cd4_ratio_pdf, 
       cd8_cd4_ratio_plot, 
       width = 2.75, height = 4, useDingbats = FALSE)

