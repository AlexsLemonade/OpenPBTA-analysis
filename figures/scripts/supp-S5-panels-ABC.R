# Stephanie J. Spielman for CCDL, 2022

# This script creates three S5 panels (A, B, C)


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
ggsave(roc_file, roc_plot, width = 5, height = 5) 


## Figures S5B and S5C---------------------------------------------------
# TERT (S5B) and TERC (S5C) expression across normextend scores

# Paths and data for this figure (and the next)
telomerase_dir <- file.path(analyses_dir, "telomerase-activity-prediction")
# stranded data specifically:
extend_scores <- read_tsv(file.path(telomerase_dir, "results", "TelomeraseScores_PTBAStranded_FPKM.txt")) 
stranded_expression <- read_rds(file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"))

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
  mutate(FPKM = log(FPKM + 1))

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
    aes(x = NormEXTENDScores, 
        y = FPKM) + 
    # All circular shapes (1, 16, 19, 20, 21) are appearing as greek letter lambdas during PDF panel compilation
    # But, squares seem to work!
    geom_point(shape = 15) + 
    geom_smooth(method = "lm") + 
    annotate("text", 
             label = stats_df$annotation[stats_df$gene == gene_name], 
             x = 0.2, 
             y = annotation_y,
             size = 3) + 
    labs(x = "Telomerase score",
         y = paste0(gene_name, " log(FPKM)")
    ) +
    ggpubr::theme_pubr()
}

tert_plot <- plot_extend_scatter(extend_fpkm_df, stats_annotation_df, "TERT", 4.5) #S5b
terc_plot <- plot_extend_scatter(extend_fpkm_df, stats_annotation_df, "TERC", 6) #S5c


ggsave(tert_file, tert_plot, width = 5, height = 5)
ggsave(terc_file, terc_plot, width = 5, height = 5)






