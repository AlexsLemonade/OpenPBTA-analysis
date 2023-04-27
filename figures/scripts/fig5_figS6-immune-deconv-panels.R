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


## Input and output files -------------------------------------

hist_file <- file.path(data_dir, "pbta-histologies.tsv")
palette_file <- file.path(palette_dir, "broad_histology_cancer_group_palette.tsv")
quantiseq_file <- file.path(analysis_dir, "quantiseq_deconv-output.rds")
polya_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")
stranded_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

# Figure 5 output:
quantiseq_cancer_group_pdf <- file.path(output_dir_fig5, "quantiseq-cell_types-cancer_groups.pdf")
cd274_expression_mb_pdf <- file.path(output_dir_fig5, "cd274_expression_mb_subtypes.pdf")

# Figure S6 output:
quantiseq_subtypes_pdf <- file.path(output_dir_figS6, "quantiseq-cell_types-molecular_subtypes.pdf")
cd8_cd4_ratio_pdf <- file.path(output_dir_figS6, "cd8_cd4_ratio.pdf")

# Zenodo CSV output directory and file paths
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig5c_csv <- file.path(zenodo_tables_dir, "figure-5c-data.csv")
fig5e_csv <- file.path(zenodo_tables_dir, "figure-5e-data.csv")
figS6e_csv <- file.path(zenodo_tables_dir, "figure-S6e-data.csv")
figS6f_csv <- file.path(zenodo_tables_dir, "figure-S6f-data.csv")



# Prepare data for plots --------------

# Read in clinical data and palette file
histologies_df <- read_tsv(hist_file, guess_max = 10000)
palette_df <- readr::read_tsv(palette_file)

# Read in analysis data
quantiseq <- read_rds(quantiseq_file) %>%
  # Unneeded columns
  select(-library, -method)
polya <- read_rds(polya_file)
stranded <- read_rds(stranded_file)

# Prepare palette/histology data
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
         molecular_subtype) %>%
  # Change subtype loss -> lost
  mutate(molecular_subtype = str_replace(molecular_subtype, "loss", "lost"))





## Plot 5C ----------------------


# Cancer groups we want to visualize in 5C
cancer_groups_of_interest_5c <- c("Pilocytic astrocytoma",
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

# filter to relevant samples and remove uncharacterized fractions
quantiseq_cg <- quantiseq %>%
  filter(cell_type != "uncharacterized cell") %>%
  left_join(palette_mapping_df, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  filter(cancer_group_display %in% cancer_groups_of_interest_5c)

# We need to create a new scheme for labeling that shows wrapped cancer groups with `(n=X)`
quantiseq_cg <- quantiseq_cg %>%
  count(cancer_group_display, cell_type) %>%
  select(-cell_type) %>%
  distinct() %>%
  inner_join(
    select(palette_df, cancer_group_display, cancer_group_hex)
  ) %>%
  select(cancer_group_display, n, cancer_group_hex) %>%
  # Create wrapped with (n=X) factor column for cancer groups
  mutate(cancer_group_display_n = stringr::str_wrap(glue::glue("{cancer_group_display} (N={n})"), 20),
                cancer_group_display_n = forcats::fct_reorder(cancer_group_display_n, n, .desc=T)) %>%
  inner_join(quantiseq_cg)

# plot by cancer group
cancer_group_plot <- ggplot(quantiseq_cg) +
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
    y = "Estimated fraction in tumor"
  ) +
  ggpubr::theme_pubr() +
  cowplot::panel_border() +
  theme(
    # Sizing for compilation
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = rel(0.5),
                              margin = margin(0.2, 0.25, 0.2, 0.25)), 
    strip.background = element_rect(size = 0.4),
    axis.line = element_line(size = rel(0.4)),
    axis.ticks = element_line(size = rel(0.4)),
    legend.position = "none"
  )

# Sized for compilation
ggsave(quantiseq_cancer_group_pdf,
       cancer_group_plot,
       width = 9, height = 3, useDingbats=FALSE)


## Plot 5E ------------------------------------


# Combine polya and stranded to get expression for CD274
polya_expression <- polya %>%
  rownames_to_column("gene") %>%
  filter(gene == "CD274") %>%
  gather(-gene, key = sample, value = expression)

expression_CD274 <- stranded %>%
  rownames_to_column("gene") %>%
  filter(gene == "CD274") %>%
  gather(-gene, key = sample, value = expression) %>%
  # combine with polya
  bind_rows(polya_expression) %>%
  as_tibble() %>%
  # log2(x+1) transformation
  mutate(expression = log(expression +1, 2))

# Prepare data for plotting:
CD274_cd8_mb <- expression_CD274 %>%
  # Join in CD8+ fractions
  inner_join(
    filter(quantiseq, cell_type == "T cell CD8+")
  ) %>%
  # Join in classified MB subtypes
  inner_join(
    palette_mapping_df %>%
      filter(cancer_group_display == "Medulloblastoma",
             !(str_detect(molecular_subtype, "To be classified"))
             ) %>%
      select(sample = Kids_First_Biospecimen_ID, molecular_subtype)
  ) %>%
  # Count group sizes to use in x-axis labels
  group_by(molecular_subtype) %>%
  mutate(subtype_count = n()) %>%
  ungroup() %>%
  mutate(molecular_subtype = glue::glue("{molecular_subtype}\n(N = {subtype_count})"))

# Prepare data for labeling P-values
wilcox_df <- ggpubr::compare_means(expression ~ molecular_subtype,
                                   data = CD274_cd8_mb,
                                   method = "wilcox.test") %>%
  # Add y-axis positions for P-values
  mutate(y.position = c(0.8,
                        NA,
                        1.8,
                        1.0,
                        2.0,
                        2.2))

# Plot 5E
cd274_expression_mb_plot <- ggplot(CD274_cd8_mb) +
  aes(x = molecular_subtype,
      # Note this is log2(fpkm+1)!
      y = expression) +
  # remove outliers
  geom_boxplot(outlier.shape = NA,
               color = "grey20",
               size = 0.2) +
  geom_jitter(width = 0.15,
              size = 0.75,
              alpha = 0.5,
              # Helpful shape for compiling
              shape = 16) +
  ggpubr::stat_pvalue_manual(wilcox_df,
                             label = "p = {p.adj}",
                             label.size = 2,
                             bracket.size = 0.25) +
  scale_color_identity() +
  labs(x = "Molecular subtype of tumor",
       y = "CD274 log2(FPKM+1)") +
  ggpubr::theme_pubr() +
  # Set sizing for compilation
  theme(
    axis.text.x = element_text(hjust = 0.55,
                               size = 6.5),
    axis.text.y = element_text(size = 6.5),
    axis.title = element_text(size = 7),
    axis.line  = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2)
    )

# Panel sized for compilation
ggsave(cd274_expression_mb_pdf,
       cd274_expression_mb_plot, width = 3, height = 2.25,
       useDingbats=FALSE)



## Plot S6E -----------------------------------


# Histologies to include in this panel and in the next panel S6F - https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/1575#discussion_r930516213
#  But ensure "Ependymoma" terminology is used: https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1652
broad_histology_order <- c("Tumor of sellar region", "Ependymoma", "Embryonal tumor", "High-grade glioma")


# Find molecular subtypes, and their order, to include in this panel
subtype_order <- palette_mapping_df %>%
  select(broad_histology_display, molecular_subtype) %>%
  # keep only relevant histologies
  filter(broad_histology_display %in% broad_histology_order) %>%
  # remove NA and unclassified subtypes
  filter(!is.na(molecular_subtype),
         !str_detect(molecular_subtype, "To be classified")) %>%
  # Keep only combinations with N>=3
  count(broad_histology_display, molecular_subtype) %>%
  filter(n >= 3) %>%
  # Factor/arrange broad_histology_display to obtain the final molecular_subtype order
  mutate(broad_histology_display = fct_relevel(broad_histology_display, broad_histology_order)) %>%
  arrange(broad_histology_display) %>%
  pull(molecular_subtype)



# Establish data for this plot
data_for_s6e <- quantiseq %>%
  # remove uncharacterized cells
  filter(cell_type != "uncharacterized cell") %>%
  # Join in palette/subtype information
  inner_join(
    select(palette_mapping_df,
           sample = Kids_First_Biospecimen_ID,
           broad_histology_display, broad_histology_hex, molecular_subtype)
  ) %>%
  # Filter and order subtypes, broad histology
  filter(molecular_subtype %in% subtype_order) %>%
  mutate(molecular_subtype = fct_relevel(molecular_subtype, subtype_order),
         broad_histology_display = fct_relevel(broad_histology_display, broad_histology_order))

# Set up palette for this plot
pal_df <- data_for_s6e %>%
  select(broad_histology_display, broad_histology_hex) %>%
  distinct() %>%
  arrange(broad_histology_display)
bh_hex <- pal_df$broad_histology_hex
names(bh_hex) <- pal_df$broad_histology_display


# Create the plot
quantiseq_subtypes_plot <- ggplot(data_for_s6e) +
  aes(x = molecular_subtype, y = score) +
  geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.2) +
  geom_jitter(width = 0.15, size = 0.66, alpha = 0.6,
              aes(color = broad_histology_display),
              # Helpful shape for compiling
              shape = 16) +
  facet_wrap(~cell_type, ncol = 5, scales = "free_y") +
  scale_color_manual(values = bh_hex, name = "Broad histology") +
  labs(
    x = "Molecular subtype of tumor",
    y = "Estimated fraction in tumor"
  ) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.title = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(0.7)),
        strip.text = element_text(size = rel(0.7)),
        axis.line = element_line(size = rel(0.7)),
        axis.ticks = element_line(size = rel(0.7)),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.7))) +
  # set a larger, visible legend point size without alpha
  guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))


ggsave(quantiseq_subtypes_pdf,
       quantiseq_subtypes_plot,
       width = 11, height = 4, useDingbats=FALSE)




## Plot S6F -----------------------------------

# Calculate cd8+/cd4+ ratio for all subtypes, join with subtype/histology, and set up factors
ratio_df <- quantiseq %>%
  spread(cell_type, score) %>%
  mutate(cd8_cd4_ratio = `T cell CD8+` / `T cell CD4+ (non-regulatory)`)  %>%
  gather(cell_type, score, -sample) %>%
  # Keep only the known ratios
  filter(cell_type == "cd8_cd4_ratio",
         !(is.infinite(score)),
         !(is.nan(score))) %>%
  rename(cd8_cd4_ratio = cell_type) %>%
  # Join in palette/subtype information
  inner_join(
    select(palette_mapping_df,
           sample = Kids_First_Biospecimen_ID,
           broad_histology_display, broad_histology_hex, molecular_subtype)
  ) %>%
  # Filter to subtype_order set up for Figure S6E
  filter(molecular_subtype %in% subtype_order)

# Update the subtype order to only keep those that are in the ratio_df
#  Some may not have been retained if they were all NaN or Inf
subtype_order <- subtype_order[subtype_order %in% unique(ratio_df$molecular_subtype)]

# Similarly update the histology order to only keep those that are in the ratio_df
broad_histology_order <- broad_histology_order[broad_histology_order %in% unique(ratio_df$broad_histology_display)]


# Update factor order for subtypes and histology and plot
cd8_cd4_ratio_plot <- ratio_df %>%
  mutate(molecular_subtype = fct_relevel(molecular_subtype, subtype_order),
         broad_histology_display = fct_relevel(broad_histology_display, broad_histology_order)) %>%
  ggplot() +
  aes(x = molecular_subtype,
      y = score) +
  # remove outliers
  geom_boxplot(outlier.shape = NA, color = "grey40", size = 0.5) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6,
              aes(color = broad_histology_hex),
              # Helpful shape for compiling
              shape = 16) +
  scale_color_identity() +
  labs(x = "Molecular subtype of tumor",
       y = "Ratio of CD8+/CD4+ T cell fractions") +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(hjust = 1,
                                   size = rel(0.6),
                                   angle = 45),
        axis.text.y = element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7)),
        axis.line = element_line(size = rel(0.7)),
        axis.ticks = element_line(size = rel(0.7)))


ggsave(cd8_cd4_ratio_pdf,
       cd8_cd4_ratio_plot,
       width = 3.5, height = 4, useDingbats = FALSE)



# Export CSVs for Zenodo upload

# Panel 5C
quantiseq_cg %>%
  # reorder columns so the ID is first
  # and remove column with \n!
  dplyr::select(Kids_First_Biospecimen_ID = sample, everything(), -cancer_group_display_n) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  # export
  readr::write_csv(fig5c_csv)


# Panel 5E
CD274_cd8_mb %>%
  # reorder columns so the ID is first
  dplyr::select(Kids_First_Biospecimen_ID = sample, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  # remove \n from molecular_subtype so that CSV is properly formatted
  dplyr::mutate(molecular_subtype = stringr::str_replace(molecular_subtype, "\n.+", "")) %>%
  # export
  readr::write_csv(fig5e_csv)





# Panel S6E
data_for_s6e %>%
  # reorder columns so the ID is first
  dplyr::select(Kids_First_Biospecimen_ID = sample, everything(),
                # but remove the hex column
                -broad_histology_hex) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  # export
  readr::write_csv(figS6e_csv)




# Panel S6F
ratio_df %>%
  # reorder columns so the ID is first
  dplyr::select(Kids_First_Biospecimen_ID = sample, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  # export
  readr::write_csv(figS6f_csv)


