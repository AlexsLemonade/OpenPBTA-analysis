# Stephanie J. Spielman for CCDL, 2022-3

# This script creates panels for Figure S5


## Load libraries ------------------------
library(tidyverse)


## Define paths  -----------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs5", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
analyses_dir <- file.path(root_dir, "analyses")
data_dir <- file.path(root_dir, "data")

# Palette file
binary_palette_df <- read_tsv(file.path(root_dir, "figures", "palettes", "binary_color_palette.tsv"))

## Define output files
roc_file <- file.path(output_dir, "supp_roc_tp53_polya.pdf") ## ROC curve
terc_file <- file.path(output_dir, "supp_terc_normextend.pdf") ## TERC vs normextend
tert_file <- file.path(output_dir, "supp_tert_normextend.pdf") ## TERT vs normextend

# Zenodo CSV output directory and file paths
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS5a_csv <- file.path(zenodo_tables_dir, "figure-S5a-data.csv")
figS5b_csv <- file.path(zenodo_tables_dir, "figure-S5b-data.csv")
figS5c_csv <- file.path(zenodo_tables_dir, "figure-S5c-data.csv")




## Figure S5A ---------------------------------------------------
# ROC curve for tp53 classifier specifically of polyA samples

# Paths and data for this figure
tp53_dir <- file.path(analyses_dir, "tp53_nf1_score")
tp53_roc_polya       <- read_tsv(file.path(tp53_dir, "results", "polya_TP53_roc_threshold_results.tsv"))
tp53_roc_polya_shuff <- read_tsv(file.path(tp53_dir, "results", "polya_TP53_roc_threshold_results_shuffled.tsv"))

# Prep ROC data
roc_df <- bind_rows(tp53_roc_polya, tp53_roc_polya_shuff) %>%
  mutate(auroc = round(auroc, 2),
         shuffled = ifelse(shuffled, 'TP53 Shuffle', 'TP53'),
         Classifier = paste0(shuffled, ' (AUROC = ', auroc,')'))

# Prep binary color palette
binary_scale <- binary_palette_df$hex_codes[binary_palette_df$color_names != "na_color"]


# Make the ROC curve
roc_plot <- ggplot(roc_df) + 
  aes(
    x = fpr,  
    y = tpr
  ) +
  geom_step(
    aes(color = Classifier), 
    size = 0.6
  ) + 
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1), 
    color = "black"
  ) +
  coord_fixed() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = binary_scale) +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate") + 
  ggpubr::theme_pubr() +
  theme(legend.text = element_text(size = rel(0.525)),
        legend.title = element_text(size = rel(0.525)),
        legend.key.size = unit(5, "points"),
        legend.margin = margin(0, 0, 0, -30),
        axis.text = element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7))
  )
ggsave(roc_file, roc_plot, width = 3, height = 3) 


## Figures S5B and S5C---------------------------------------------------
# TERT (S5B) and TERC (S5C) expression across normextend scores with TERTp-positive samples emphasized

# Paths and data for this figure (and the next)
telomerase_dir <- file.path(analyses_dir, "telomerase-activity-prediction")

# stranded data specifically:
extend_scores <- read_tsv(file.path(telomerase_dir, "results", "TelomeraseScores_PTBAStranded_FPKM.txt")) 
stranded_expression <- read_rds(file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"))

# MAF files to find TERTp-positive samples
maf_df <- read_tsv(file.path(data_dir, "pbta-snv-consensus-mutation.maf.tsv.gz"))
hotspot_df <- read_tsv(file.path(data_dir, "pbta-snv-scavenged-hotspots.maf.tsv.gz"))

# Metadata
metadata_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# Determine which samples have TERTp
snps <- c("rs1561215364", "rs1242535815")
maf_snp_bio_ids <- maf_df %>%
  filter(dbSNP_RS %in% snps) %>%
  pull(Tumor_Sample_Barcode) 

hotspot_snp_bio_ids <- hotspot_df %>%
  filter(dbSNP_RS %in% snps) %>%
  pull(Tumor_Sample_Barcode)

# There are 9 samples:
snp_biospecimen_ids <- unique(c(maf_snp_bio_ids, hotspot_snp_bio_ids))

# Get the corresponding sample ids
tertp_sample_ids <- metadata_df %>%
  filter(Kids_First_Biospecimen_ID %in% snp_biospecimen_ids) %>%
  pull(sample_id) 

# Finally, get the corresponding biospecimen IDs for those sample ids 
#  which have stranded RNA libraries
rna_tertp_biospecimen_ids <- metadata_df %>%
  filter(sample_id %in% tertp_sample_ids, 
         experimental_strategy == "RNA-Seq",
         RNA_library == "stranded") %>%
  pull(Kids_First_Biospecimen_ID) 


# Combine extend scores with expression for genes of interest
extend_fpkm_df <- stranded_expression %>% 
  rownames_to_column("gene") %>% 
  # keep only genes we are visualizing
  filter(gene %in% c("TERC", "TERT")) %>%
  gather(contains("BS"), key = "SampleID", value = "FPKM") %>%
  inner_join(
    select(extend_scores, 
           SampleID,
           NormEXTENDScores), 
    by = "SampleID"
  ) %>%
  mutate(FPKM = log(FPKM + 1, 2),
         # Indicate TERTp samples
         tertp = ifelse(SampleID %in% rna_tertp_biospecimen_ids, 
                        "TERTp mutation present", 
                        "TERTp mutation not observed")
  )

# Calculate stats
extend_fpkm_lm <- function(df) {
  lm(FPKM ~ NormEXTENDScores, data = df)
}

stats_annotation_df <- extend_fpkm_df %>%
  # build a model per gene
  group_by(gene) %>%
  nest(FPKM, NormEXTENDScores) %>%
  mutate(
    # build lms and then glance to get R^2 and P-value
    fit = map(data, extend_fpkm_lm),
    glanced = map(fit, broom::glance)
  ) %>%
  select(gene, glanced) %>%
  # Reveal glanced output
  unnest(glanced) %>%
  # select stats columns of interest
  select(gene, r.squared, p.value) %>%
  # Format a column to use for figure annotation
  mutate(r = sqrt(r.squared),
         annotation = paste0("R = ", 
                             round(r,3), 
                             "; P-value = ", 
                             format(p.value, digits=3)))

plot_extend_scatter <- function(plot_df, stats_df, gene_name, annotation_y) {
  plot_df %>%
    filter(gene == gene_name) %>%
    ggplot() + 
    geom_point(aes(color = tertp, 
                   alpha = tertp, 
                   size = tertp)) +
    aes(x = NormEXTENDScores, 
        y = FPKM) + 
    geom_smooth(method = "lm", 
                color = "black") + 
    scale_color_manual(values = c("grey70", "red")) +
    scale_alpha_manual(values = c(0.5, 1), guide = FALSE) +
    scale_size_manual(values = c(rel(1), rel(1.3)), guide = FALSE) +
    annotate("text", 
             label = stats_df$annotation[stats_df$gene == gene_name], 
             x = 0.28, 
             y = annotation_y,
             size = 2.5) + 
    labs(x = "Telomerase score",
         y = paste0(gene_name, " log2(FPKM+1)"),
         color = ""
    ) +
    ggpubr::theme_pubr() + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6))
}

tert_plot <- plot_extend_scatter(extend_fpkm_df, stats_annotation_df, "TERT", 6.5) 
terc_plot <- plot_extend_scatter(extend_fpkm_df, stats_annotation_df, "TERC", 8.5) 


ggsave(tert_file, tert_plot, width = 4, height = 4, useDingbats = FALSE)
ggsave(terc_file, terc_plot, width = 4, height = 4, useDingbats = FALSE)



## Export CSVs for Zenodo upload ------------------------------

# Panel S5A: ROC curve
# no sample information so no arranging is needed
readr::write_csv(roc_df, figS5a_csv)


# Panel S5B and S5C are made from the same data, 
#  so prep first and then filter for each gene
extend_fpkm_df_export <- extend_fpkm_df %>%
  dplyr::select(Kids_First_Biospecimen_ID = SampleID, everything()) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID)


# S5B
extend_fpkm_df_export %>%
  dplyr::filter(gene == "TERT") %>%
  readr::write_csv(figS5b_csv)


# S5C
extend_fpkm_df_export %>%
  dplyr::filter(gene == "TERC") %>%
  readr::write_csv(figS5c_csv)





