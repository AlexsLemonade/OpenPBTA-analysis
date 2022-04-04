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

# Read in TMB for highlighting points in boxplots
tmb_coding_df <- read_tsv(file.path(data_dir, "pbta-snv-consensus-mutation-tmb-coding.tsv"))
                            

#### Define output PDF panels ----------------------------------------------------------------
tp53_roc_pdf                        <- file.path(output_dir, "tp53_stranded_roc_panel.pdf")
tp53_scores_altered_pdf             <- file.path(output_dir, "tp53_scores_by_altered_panel.pdf") 
tp53_expression_altered_pdf         <- file.path(output_dir, "tp53_expression_by_altered_panel.pdf") 
tp53_telomerase_scores_boxplot_pdf        <- file.path(output_dir, "tp53_telomerase_boxplots_panel.pdf") 
tp53_telomerase_scores_boxplot_legend_pdf <- file.path(output_dir, "tp53_telomerase_boxplots_panel_legend.pdf")
#forest_plot_pdf                     <- file.path(output_dir, "forest_survival_tp53_telomerase_hgg_panel.pdf") 



### ROC curve ---------------------------------------------------------------------------------

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




 ### TP53 scores and expression violin plots -----------------------------------------------------------
# We do not use color palettes since the color mappings will necessarily change across figures.
# For all violin plots, we will just show `activated` and `loss` because `other` is 
#    actually just **unclassified** and therefore not a robust comparison group
# We use ggplot instead of ggpubr because of jitter styling (ggpubr jitter point placement is deterministic)


# remove "other" grouping from data so we only have 2 categories and 
#  change to `activated` and `lost` for grammatical consistency
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
ggsave(tp53_scores_altered_pdf, tp53_scores_plot, width = 6, height = 4)
ggsave(tp53_expression_altered_pdf, tp53_expression_plot, width = 6, height = 4)




## TP53 and telomerase scores boxplots across cancer groups with mutators emphasized -------------------------------------------


# Define cancer groups to show in boxplots
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

# Cancer group wrap number of characters for labeling x-axis in boxplots
cg_wrap <- 20


# Get coding samples of interest and annotate samples as normal, hyper, or ultra
mutators <- tmb_coding_df %>%
  # columns of interest
  select(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode,
         tmb) %>%
  # keep only rows of interest
  filter(Kids_First_Biospecimen_ID %in% tp53_scores_df_subset$Kids_First_Biospecimen_ID) %>%
  # add a column about mutator
  mutate(
    mutator = case_when(
      tmb < 10 ~ "Normal",
      tmb >= 10 & tmb < 100 ~ "Hypermutant",
      tmb >= 100 ~ "Ultra-hypermutant")
  ) %>%
  # find RNA id version
  inner_join(
    select(histologies_df, 
             Kids_First_Biospecimen_ID)
  )

# Combine extend (telomerase) and tp53 scores
extend_tp53_df <- extend_scores %>%
  # consistenct naming for joining
  select(Kids_First_Biospecimen_ID = SampleID,
         NormEXTENDScores) %>%
  # join with tp53
  left_join(
    select(tp53_compare,
           # consistenct naming for joining
           Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA,
           tp53_score),
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  # join with metadata
  inner_join(
    select(histologies_df, 
           Kids_First_Biospecimen_ID,
           cancer_group)
  ) %>%
  # join with palette
  inner_join(
    select(histologies_palette_df,
           cancer_group, 
           cancer_group_display,
           cancer_group_hex)
  ) %>%
  distinct() %>%
  # filter to cancer groups of interest
  filter(cancer_group %in% cancer_groups_to_plot) %>%
  # join with mutator information
  inner_join(mutators)
  
  
# Prepare combined data for visualization
extend_tp53_plot_df <- extend_tp53_df %>%
  # duplicate the tp53 scores column so we can eventually order can groups by it
  mutate(tp53_score_forordering = tp53_score) %>%
  # we want a single column for all scores so we can facet by scores
  gather(NormEXTENDScores, tp53_score, 
         key = "score_type", 
         value = "score")
  # wrap cancer group label
  mutate(cancer_group_display = str_wrap(cancer_group_display, cg_wrap)) %>%
  
  


# Join data for visualization
tp53_plot_df <- tp53_scores_df_subset %>%
  inner_join(mutators, by="Kids_First_Biospecimen_ID") %>%
  inner_join(
    select(histologies_df, 
           Kids_First_Biospecimen_ID,
           cancer_group)
  ) %>%
  # filter to cancer groups and wrap for x-axis
  filter(cancer_group %in% cancer_groups_to_plot) %>%
  mutate(cancer_group = str_wrap(cancer_group, cg_wrap))


## Make plot ans associated legend
# Define colors to use
legend_colors <- c(Normal      = "grey40",
                   Hypermutant = "orange",
                   `Ultra-hypermutant` = "red")

# Boxplot with overlayed jitter with colors "mapped" to `mutator`
set.seed(9) # reproducible jitter to ensure we can see all N=6 points
tp53_tmb_boxplot <- ggplot(tp53_plot_df) +
  aes(
    x = fct_reorder(cancer_group, tp53_score), # order cancer groups by tp53 score
    y = tp53_score
  ) +
  geom_boxplot(
    outlier.shape = NA, # no outliers
    color = "grey20",    # dark grey color
    size = 0.4 
  ) +
  # Separate out jitters so that the mutant layers are ON TOP OF normal
  geom_jitter(
    data = tp53_plot_df[tp53_plot_df$mutator == "Normal",],
    width = 0.15, 
    alpha = 0.7,
    pch = 19,
    color = legend_colors["Normal"]
  ) +
  geom_jitter(
    data = tp53_plot_df[tp53_plot_df$mutator == "Hypermutant",],
    width = 0.15, 
    pch = 21,
    size = 2.25,
    fill = legend_colors["Hypermutant"]
  ) +  
  geom_jitter(
    data = tp53_plot_df[tp53_plot_df$mutator == "Ultrahypermutant",],
    width = 0.15, 
    pch = 21,
    size = 2.25,
    fill = legend_colors["Ultrahypermutant"]
  ) +  
  labs(x = "Cancer group", 
       y = "TP53 Score") +
  # do we want an hline at 0.5? Might be useful guiding line but also add unnecessary visual noise
  # geom_hline(yintercept = 0.5) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = rel(0.8))
  )

# Export plot
ggsave(tp53_scores_cancergroups_pdf,
       tp53_tmb_boxplot,
       width = 8, height = 5)

# Make a legend for the grey/orange/red since this was not done with normal mapping
# We have to make a "fake" plot for this to extraxt the legend from
tp53_plot_legend_df <- tp53_plot_df %>%
  mutate(mutator_factor = factor(mutator, levels=names(legend_colors)))

tp53_plot_for_legend <- ggplot(tp53_plot_legend_df) + 
  aes(x = cancer_group, y = tp53_score, shape = mutator_factor, fill = mutator_factor, color = mutator_factor) +
  geom_point(size =3) + 
  scale_shape_manual(name = "Mutator", values = c(19, 21, 21)) +
  scale_color_manual(name = "Mutator",values = c(unname(legend_colors["Normal"]), "black", "black")) + # for reasons (?) this apparently needs unname(). weird since fill doesnt
  scale_fill_manual(name = "Mutator", values = c("black", legend_colors["Hypermutant"], legend_colors["Ultrahypermutant"])) +
  # theme to remove gray background. this strategy works
  theme_classic()


legend <- cowplot::get_legend(tp53_plot_for_legend)

# Export legend
pdf(tp53_scores_cancergroups_legend_pdf, width = 6, height = 3)
cowplot::ggdraw(legend)
dev.off()


## Distributions of Extend scores ------------------------------------------------------


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
  # join with palette
  inner_join(
    select(histologies_palette_df,
           contains("cancer_group"))
  ) %>%
  # filter to cancer groups of interest
  filter(cancer_group %in% cancer_groups_to_plot) %>%
  # wrap label
  mutate(cancer_group_display = str_wrap(cancer_group_display, cg_wrap))

extend_boxplot <- ggplot(extend_df) +
  aes(x = fct_reorder(cancer_group_display, NormEXTENDScores),
      y = NormEXTENDScores, 
      color = cancer_group_hex) + 
  geom_boxplot() + 
  geom_boxplot(
    outlier.shape = NA, # no outliers
    color = "grey20",    # dark grey color
    size = 0.4 
  ) +
  geom_jitter(alpha = 0.4, 
              size = 1) + 
  labs(
    x = "Cancer group",
    y = "Telomerase scores"
  ) +
  scale_color_identity() +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = rel(0.8))
  )


# Export plot
ggsave(telomerase_scores_cancergroups_pdf,
       extend_boxplot,
       width = 8, height = 5)


