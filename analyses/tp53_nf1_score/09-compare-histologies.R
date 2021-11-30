suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
palette_dir <- file.path(root_dir, "figures", "palettes")
output_dir <- file.path(root_dir, "analyses", "tp53_nf1_score",
                        "plots")

# read in histologies and tp53 scores file
hist_file <- file.path(root_dir, "data", "pbta-histologies.tsv")
tp53_file <- file.path(root_dir, "analyses", "tp53_nf1_score",
                 "results", "tp53_altered_status.tsv")

# Declare output filenames
tp53_broad_pdf <- file.path(output_dir, "tp53_broad_histology.pdf")
tp53_cg_pdf <- file.path(output_dir, "tp53_cancer_group.pdf")


# read in histologies file
histologies <- read_tsv(hist_file, guess_max = 10000) %>%
  filter(sample_type == "Tumor") %>%
  # remove cell lines from plots
  filter(composition != "Derived Cell Line") %>%
  # select only columns of interest
  select(Kids_First_Biospecimen_ID, sample_id, broad_histology, cancer_group, molecular_subtype)

# read in tp53 scores
tp53 <- readr::read_tsv(tp53_file) %>%
  mutate(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA) %>%
  select(Kids_First_Biospecimen_ID,
         tp53_score) %>%
  filter(!is.na(Kids_First_Biospecimen_ID)) %>%
  distinct()

# join histology data with tp53 scores
tp53_hist_for_broad <- tp53 %>%
  left_join(histologies) %>%
  distinct() %>%
  # remove broad histology == NA
  filter(!is.na(broad_histology)) %>%
  # remove precancerous lesion, non-tumor, benign
  filter(broad_histology != "Pre-cancerous lesion" & broad_histology != "Non-tumor" &
           broad_histology != "Other tumor")

# remove cancer group == NA
tp53_hist_for_cancer_group <- tp53 %>%
  left_join(histologies) %>%
  distinct() %>%
  filter(!is.na(cancer_group))

# color palette
color_palette_df <- readr::read_tsv(
  file.path(palette_dir, "broad_histology_cancer_group_palette.tsv")
)

# get palette for cancer group
cancer_group_palette <- color_palette_df %>%
  dplyr::select(cancer_group, cancer_group_hex) %>%
  # Remove NA values -- a cancer group hex value will be NA only if the
  # cancer group is NA
  dplyr::filter(complete.cases(.))

# make cancer group color palette suitable for use with ggplot
col_cg <- cancer_group_palette$cancer_group_hex
names(col_cg) <- cancer_group_palette$cancer_group

# repeat the steps for the broad histology palette
broad_histology_palette <- color_palette_df %>%
  dplyr::select(broad_histology, broad_histology_hex) %>%
  dplyr::distinct() %>%
  dplyr::filter(complete.cases(.))
col_broad <- broad_histology_palette$broad_histology_hex
names(col_broad) <- broad_histology_palette$broad_histology

## boxplot of TP53 scores by broad histology
broad_histology_plot <- ggplot(tp53_hist_for_broad,
                               aes(reorder(broad_histology, -tp53_score, FUN = median, na.rm = TRUE),
                                   tp53_score)) +
  geom_boxplot(aes(color = broad_histology, fill = broad_histology), alpha = 0.3) +
  geom_hline(aes(yintercept=0.5), linetype = "dashed", color = "black") +
  geom_jitter(aes(color = broad_histology, alpha = 0.2)) +
  scale_fill_manual(values = col_broad, aesthetics = c("colour", "fill")) +
  xlab("Broad histology") +
  ylab("TP53 scores") +
  theme_pubr() +
  rremove("legend") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, size = 6, vjust = 1, hjust = 1))
ggsave(filename = tp53_broad_pdf, plot = broad_histology_plot)

## boxplot of TP53 scores by cancer group
cancer_group_plot <- ggplot(tp53_hist_for_cancer_group,
                            aes(reorder(cancer_group, -tp53_score, FUN = median, na.rm = TRUE),
                                tp53_score)) +
  geom_boxplot(aes(color = cancer_group, fill = cancer_group), alpha = 0.3) +
  geom_hline(aes(yintercept=0.5), linetype = "dashed", color = "gray60") +
  geom_jitter(aes(color = cancer_group, alpha = 0.2)) +
  scale_fill_manual(values = col_cg, aesthetics = c("colour", "fill")) +
  xlab("Cancer group") +
  ylab("TP53 scores") +
  theme_pubr() +
  rremove("legend") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, size = 6, vjust = 1, hjust = 1))
ggsave(filename = tp53_cg_pdf, plot = cancer_group_plot)


