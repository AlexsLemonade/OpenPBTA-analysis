# S. Spielman for ALSF CCDL 2022
#
# Makes pdf panels for reporting TP53 and telomerase results in main text

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
tp53_dir <- file.path(analyses_dir, "tp53_nf1_score")
telomerase_dir <- file.path(analyses_dir, "telomerase-activity-prediction")

# Palette directory
palette_dir <- file.path(root_dir, "figures", "palettes")

# Read in clinical data and associated palette
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)
histologies_palette_df <- read_tsv(file.path(palette_dir, "broad_histology_cancer_group_palette.tsv"))
binary_palette_df <- readr::read_tsv(file.path(palette_dir, "binary_color_palette.tsv"))

# Read in tp53 data for ROC
tp53_roc_stranded       <- read_tsv(file.path(tp53_dir, "results", "stranded_TP53_roc_threshold_results.tsv"))
tp53_roc_stranded_shuff <- read_tsv(file.path(tp53_dir, "results", "stranded_TP53_roc_threshold_results_shuffled.tsv"))


# Read in tp53 data file for violin plots
tp53_compare <- read_tsv(file.path(tp53_dir, "results", "tp53_altered_status.tsv"))
stranded_expression <- read_rds(file.path(data_dir,"pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"))


# Read in EXTEND scores. We read Stranded FPKM here. 
extend_scores <- read_tsv(file.path(telomerase_dir, "results", "TelomeraseScores_PTBAStranded_FPKM.txt"))


#### Define output PDF panels ----------------------------------------------------------------
tp53_roc_pdf <- file.path(output_dir, "tp53_roc_panel.pdf")
tp53_scores_pdf <- file.path(output_dir, "tp53_scores_by_altered_panel.pdf") 
tp53_expression_pdf <- file.path(output_dir, "tp53_expression_by_altered_panel.pdf") 


### ROC curve ----------------------------------------------

# Create data frame that will plot ROC
roc_df <- bind_rows(tp53_roc_stranded, tp53_roc_stranded_shuff) %>%
  mutate(auroc = round(auroc, 2),
         shuffled = ifelse(shuffled, 'TP53 Shuffle', 'TP53'),
         Classifier = paste0(shuffled, ' (AUROC = ', auroc,')'))

# prep in binary color palette
binary_scale <- binary_palette_df$hex_codes[binary_palette_df$color_names != "na_color"]


# Make the ROC curve
roc_plot <- ggplot(roc_df) + 
  aes(
    x = fpr,  
    y = tpr
  ) +
  geom_step(
    aes(color = Classifier), 
    size = 0.75
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
  theme(legend.text = element_text(size = rel(0.7)),
        legend.title = element_text(size = rel(0.7)),
        axis.text = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(0.75))
  )
ggsave(tp53_roc_pdf, roc_plot, width = 5, height = 5) 




 ### TP53 scores and expression violin plots ---------------------------------
# We do not use color palettes since the color mappings will necessarily change across figures.
# For all violin plots, we will just show `activated` and `loss` because `other` is 
#    actually just **unclassified** and therefore not a robust comparison group
# We use ggplot instead of ggpubr because of jitter styling (ggpubr jitter point placement is deterministic)


# remove "other" grouping from data so we only have 2 categories and change to `activated` and `lost` for grammatical consistency
 tp53_compare_2cat <- tp53_compare %>%
  filter(tp53_altered != "other") %>%
  # change loss --> lost
  mutate(tp53_altered = ifelse(
    tp53_altered == "loss", "lost", "activated"
  ))

# Perform test for TP53 scores
scores_pvalue <- round(
  wilcox.test(tp53_score ~ tp53_altered, data = tp53_compare_2cat)$p.value,
  2)

# Plot tp53 scores with P-value annotation
set.seed(4)
tp53_scores_plot <- ggplot(tp53_compare_2cat) +
  aes(x = tp53_altered, 
      y = tp53_score) +
  geom_violin() + 
  geom_jitter(alpha = 0.5, 
              width = 0.2, 
              size = 1) + 
  # Add mean with stat_summary
  stat_summary(color = "firebrick") + 
  # Add p-value annotation
  annotate("text", 
           label = paste0("Wilcoxon P-value = ", scores_pvalue), 
           x = 1, y = 1.2) +
  labs(
    x = "TP53 alteration status",
    y = "TP53 score"
  ) +
  ggpubr::theme_pubr() 




# subset to TP53
subset_stranded <- t(stranded_expression)[,"TP53"]

# Join STRANDED expression with tp53 alteration
# Note that because this is stranded only, it has fewer data points.
# TODO: Should the tp53 score figure also be just stranded? The ROC figure is also only just stranded.
stranded_tp53 <- as.data.frame(subset_stranded) %>%
  rename(tp53_expression=subset_stranded) %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID_RNA") %>%
  # easier to work with
  as_tibble() %>%
  inner_join(tp53_compare_2cat, by = "Kids_First_Biospecimen_ID_RNA") %>%
  # keep only columns we need
  select(Kids_First_Biospecimen_ID_RNA, tp53_expression, tp53_altered) %>%
  distinct()


# Perform test for TP53 expression
expression_pvalue <- round(
  wilcox.test(tp53_expression ~ tp53_altered, data = stranded_tp53)$p.value,
  4)

# Plot tp53 expression with P-value annotation
# TODO: analyze with log'd y instead?
set.seed(4)
tp53_expression_plot <- ggplot(stranded_tp53) +
  aes(x = tp53_altered, 
      y = tp53_expression) +
  geom_violin() + 
  geom_jitter(alpha = 0.5, 
              width = 0.2, 
              size = 1.25) + 
  # Add mean with stat_summary
  stat_summary(color = "firebrick") + 
  # Add p-value annotation
  annotate("text", 
           label = paste0("Wilcoxon P-value = ", expression_pvalue), 
           x = 1, y = 75) +
  labs(
    x = "TP53 alteration status",
    y = "TP53 expression (FPKM)"
  ) +
  ggpubr::theme_pubr() 


# Export figures
ggsave(tp53_scores_pdf, tp53_scores_plot, width = 6, height = 4)
ggsave(tp53_expression_pdf, tp53_expression_plot, width = 6, height = 4)


## Distributions of Extend scores --------------------------------

# Define cancer groups to show in this plot
cancer_groups_to_plot <- c(
  "Diffuse midline glioma",
  "Low-grade glioma astrocytoma",
  "Craniopharyngioma",
  "High-grade glioma astrocytoma",
  "Ganglioglioma",
  "Medulloblastoma",
  "Meningioma",
  "Ependymoma",
  "Schwannoma",
  "Dysembryoplastic neuroepithelial tumor"
)

extend_df <- extend_scores %>%
  # rename column and keep only Norm
  select(Kids_First_Biospecimen_ID = SampleID, 
         NormEXTENDScores) %>%
  # join with metadata
  inner_join(
    select(histologies_df, 
           Kids_First_Biospecimen_ID,
           cancer_group)
  ) %>%
  # filter to cancer groups of interest
  filter(cancer_group %in% cancer_groups_to_plot)

ggplot(extend_df) +
  aes(x = fct_reorder(cancer_group, NormEXTENDScores),
      y = NormEXTENDScores) + 
  geom_boxplot() + 
  geom_jitter()

tp53_extend %>%
  gather(c(tp53_score, NormEXTENDScores), key = score_type, value = score) %>%
  mutate(
    score_type = if_else(str_detect(score_type, "tp53"), "TP53 Score", "Normalized EXTEND Score")
  ) %>%
  ggplot() + 
  aes(x = fct_reorder(broad_histology_display, broad_histology_order), 
      y = score, 
      fill = broad_histology_hex) + 
  geom_violin(scale = "width", 
              alpha = 0.4) + 
  geom_jitter(width = 0.08, size = 0.15) +
  facet_wrap(~score_type, nrow=2) +
  labs(x = "Broad histology group",
       y = "") +
  scale_fill_identity() +
  ggpubr::theme_pubr() + 
  theme(
    axis.text.x = element_text(angle = 90, size = 8, hjust=1)
  )












