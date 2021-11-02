suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# set up directories
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 
root_dir <- "~/OpenPBTA-analysis/"
palette_dir <- file.path(root_dir, "figures", "palettes")
output_dir <- file.path(root_dir, "analyses", "tp53_nf1_score",
                        "plots")

# read in histologies and tp53 scores file
hist_file = file.path(root_dir, "data", "pbta-histologies.tsv")
tp53_file = file.path(root_dir, "analyses", "tp53_nf1_score",
                 "results", "tp53_altered_status.tsv")

# Declare output filenames
tp53_broad_png <- file.path(output_dir, "tp53_broad_histology.png")
tp53_cg_png <- file.path(output_dir, "tp53_cancer_group.png")


# Import standard color palettes for project
histology_label_mapping <- read_tsv(file.path(palette_dir, "histology_label_color_table.tsv"))

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


# Join histology and color mappings and labels
colors <- histologies %>%
  inner_join(histology_label_mapping, by = 
                      c("Kids_First_Biospecimen_ID",
                        "cancer_group", "broad_histology")) %>%
  filter(!is.na(cancer_group) & !is.na(broad_histology)) %>%
  select(cancer_group, cancer_group_hex_codes, broad_histology, hex_codes) %>%
  distinct()
  
# Make color key specific to these samples
col_cg <- colors$cancer_group_hex_codes
names(col_cg) <- colors$cancer_group

col_broad <- colors$hex_codes
names(col_broad) <- colors$broad_histology



## boxplot of TP53 scores by broad histology
png(tp53_broad_png, width = 8, height = 5, units = "in", res = 1200)
ggplot(tp53_hist_for_broad, aes(reorder(broad_histology, -tp53_score, FUN = median, na.rm = TRUE), tp53_score)) +
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
dev.off()

## boxplot of TP53 scores by cancer group
png(tp53_cg_png, width = 8, height = 5, units = "in", res = 1200)
ggplot(tp53_hist_for_cancer_group,  aes(reorder(cancer_group, -tp53_score, FUN = median, na.rm = TRUE), tp53_score)) +
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
dev.off()


