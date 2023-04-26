# S. Spielman for ALSF CCDL, Jo Lynne Rokita for D3b 2022
#
# Makes pdf panels for reporting TP53 and telomerase results in main text

library(tidyverse)
library(survival) # needed to parse model RDS

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
survival_dir <- file.path(analyses_dir, "survival-analysis", "results", "tp53_telomerase")

# Palette directory
palette_dir <- file.path(root_dir, "figures", "palettes")


# Zenodo CSV output directory and file paths
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig4a_csv <- file.path(zenodo_tables_dir, "figure-4a-data.csv")
fig4b_csv <- file.path(zenodo_tables_dir, "figure-4b-data.csv")
fig4c_csv <- file.path(zenodo_tables_dir, "figure-4c-data.csv")
fig4d_csv <- file.path(zenodo_tables_dir, "figure-4d-data.csv")
fig4f_csv <- file.path(zenodo_tables_dir, "figure-4f-data.csv")



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

# Read in survival model results for forest plot. 
#  Note that this object is a list:
#   The model itself is in `$model`
#   The data that went into the model is in `$data` 
survival_result <- readRDS(file.path(
  survival_dir,
  "cox_additive_terms_tp53_telomerase_resect_glioma_group.RDS"
))

#### Define output PDF panels ----------------------------------------------------------------
tp53_roc_pdf                        <- file.path(output_dir, "tp53_stranded_roc_panel.pdf")
tp53_scores_altered_pdf             <- file.path(output_dir, "tp53_scores_by_altered_panel.pdf")
tp53_expression_altered_pdf         <- file.path(output_dir, "tp53_expression_by_altered_panel.pdf")
tp53_telomerase_scores_boxplot_pdf        <- file.path(output_dir, "tp53_telomerase_boxplots_panel.pdf")
tp53_telomerase_scores_boxplot_legend_pdf <- file.path(output_dir, "tp53_telomerase_boxplots_panel_legend.pdf")
survival_plot_pdf                     <- file.path(output_dir, "forest_survival_tp53_telomerase_panel.pdf")


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
    size = 0.5
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1),
    color = "black",
    size = 0.25,
  ) +
  coord_fixed() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = binary_scale) +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate") +
  ggpubr::theme_pubr() +
  guides(color = guide_legend(title.position = "top",
                              nrow = 2)) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.title = element_text(size = rel(0.5)),
        legend.text = element_text(size = rel(0.48)),
        legend.title = element_text(size = rel(0.525)),
        legend.key.size = unit(8, "points"),
        axis.line = element_line(size = rel(0.4)),
        axis.ticks = element_line(size = rel(0.4))
  )
ggsave(tp53_roc_pdf, roc_plot, width = 1.75, height = 2,
       # add for figure compilation
       useDingbats = FALSE)



### TP53 scores and expression violin plots -----------------------------------------------------------
# We do not use color palettes since the color mappings will necessarily change across figures.
# We use ggplot instead of ggpubr because of jitter styling (ggpubr jitter point placement is deterministic)


# Function to plot both TP53 violin plots (score across altered and expression across altered)
plot_tp53 <- function(df, pvalue_y) {
  # df: Assumes two columns: the variable of interest (either tp53_score, tp53_expression) is named tp53, and column tp53_altered with three values
  # pvalue_y: y-axis coordinate of p-value annotation. Note that the x-axis coordinate is always the same
  # Return both the plot as well as data frame used to make the plot to save for Zenodo CSVs

  # seed for jitter
  set.seed(4)

  # Count group sizes to use in x-axis labels. Use this data frame going forward
  df_counts <- df %>%
    group_by(tp53_altered) %>%
    mutate(altered_counts = n()) %>%
    ungroup() %>%
    mutate(tp53_altered = glue::glue("{tp53_altered}\n(N = {altered_counts})"))


  # Perform test for variable without `other`
  df_2cat <- filter(df_counts,
                    !(str_detect(tp53_altered,"other")))

  # Prepare for use with `stat_pvalue_manual()`
  wilcox_df <- ggpubr::compare_means(tp53 ~ tp53_altered,
                                     data = df_2cat,
                                     method = "wilcox.test") %>%
    mutate(y.position = pvalue_y)



  # Prepare stats df for median +/ IQR
  stats_df <- df_counts %>%
    group_by(tp53_altered) %>%
    summarize(
      y = median(tp53, na.rm=TRUE),
      ymin = quantile(tp53, 0.25, na.rm=TRUE),
      ymax = quantile(tp53, 0.75, na.rm=TRUE)
    )

  df_counts_plot <- ggplot(df_counts) +
    aes(x = tp53_altered,
        y = tp53) +
    geom_violin(size = 0.25) +
    geom_jitter(alpha = 0.25, # very light alpha to accomodate `other` category
                width = 0.1,
                size = 0.45) +
    # Add median +/- IQR pointrange
    geom_pointrange(data = stats_df,
                    aes(
                      x = tp53_altered,
                      y = y,
                      ymin = ymin,
                      ymax = ymax
                    ),
                    color = "firebrick", size = rel(0.2)
    ) +
    # Add p-value annotation with ggpubr
    ggpubr::stat_pvalue_manual(
      wilcox_df,
      label = "Wilcoxon P-value = {p.adj}",
      # as needed, further adjust text size during compilation 
      # since larger than this ends up outside plot margins
      size = 1.75, 
      bracket.size = 0.2
    ) +
    ggpubr::theme_pubr()  +
    # Sizing for compilation - small figure export
    theme(
      axis.text = element_text(size = rel(0.7)),
      axis.title = element_text(size = rel(0.7)),
      axis.line = element_line(size = rel(0.4)),
      axis.ticks = element_line(size = rel(0.4))
    )

  return(
    list(
      "plot" = df_counts_plot,
      "df" = df_counts
    )
  )

}

# change `loss` --> `lost` for grammatical consistency
tp53_compare <- tp53_compare %>%
  mutate(tp53_altered = case_when(
    tp53_altered == "loss" ~ "lost",
    TRUE ~ tp53_altered
  ))

# Prepare data for all plots - we want stranded ONLY
# subset to TP53
subset_stranded <- t(stranded_expression)[,"TP53"]
# Join STRANDED expression with tp53 alteration
# Note that because this is stranded only, it has fewer data points.
stranded_tp53 <- as.data.frame(subset_stranded) %>%
  rename(tp53_expression=subset_stranded) %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID_RNA") %>%
  # easier to work with
  as_tibble() %>%
  # inner_join ensures stranded-only
  inner_join(tp53_compare, by = "Kids_First_Biospecimen_ID_RNA") %>%
  # keep only columns we need
  select(Kids_First_Biospecimen_ID_RNA, tp53_expression, tp53_altered, tp53_score) %>%
  distinct()


# Make the figures and obtain the tables for Zenodo CSVs
tp53_scores_plot_data <- plot_tp53(
  rename(stranded_tp53, tp53 = tp53_score),
  pvalue_y = 1.05
)
tp53_scores_data <- tp53_scores_plot_data$df
tp53_scores_plot <- tp53_scores_plot_data$plot +
  # add labels for this plot
  labs(
    x = "TP53 altered status",
    y = "TP53 score"
  )

tp53_expression_plot_data <- plot_tp53(
  stranded_tp53 %>%
    rename(tp53 = tp53_expression) %>%
    mutate(tp53 = log(tp53 + 1)),
  pvalue_y = 4.5
)
tp53_expression_data <- tp53_expression_plot_data$df

tp53_expression_plot <- tp53_expression_plot_data$plot +
  # add labels for this plot
  labs(
    x = "TP53 altered status",
    y = "TP53 expression [log(FPKM)]"
  )


# Export figures, with `useDingbats = FALSE` needed for compiling panels in Illustrator
ggsave(tp53_scores_altered_pdf, tp53_scores_plot, width = 3, height = 2, useDingbats = FALSE)
ggsave(tp53_expression_altered_pdf, tp53_expression_plot, width = 3, height = 2, useDingbats = FALSE)


## TP53 and telomerase scores boxplots across cancer groups with mutators emphasized -------------------------------------------


# Define cancer groups to show in boxplots
cancer_groups_to_plot <- c("Diffuse midline glioma",
                        "Other high-grade glioma",
                        "Pilocytic astrocytoma",
                        "Ganglioglioma",
                        "Pleomorphic xanthoastrocytoma",
                        "Other low-grade glioma",
                        "Medulloblastoma",
                        "Atypical Teratoid Rhabdoid Tumor",
                        "Other embryonal tumor",
                        "Ependymoma",
                        "Craniopharyngioma",
                        "Meningioma")

# Cancer group wrap number of characters for labeling x-axis in boxplots
cg_wrap <- 20

# Create combined data frame of mutator information, tp53 scores, and telomerase scores (NormEXTENDScores)
# Join tmb_coding_df with tp53_compare first, because both have DNA identifiers.
# Since tp53_compare also has the RNA identifier,then we can join with extend_scores
tp53_telo_mutator_df <- tmb_coding_df %>%
  # columns of interest
  select(Kids_First_Biospecimen_ID_DNA = Tumor_Sample_Barcode, # consistent naming for joining
         tmb) %>%
  # add a column about mutator
  mutate(
    mutator = case_when(
      tmb < 10 ~ "Normal",
      tmb >= 10 & tmb < 100 ~ "Hypermutant",
      tmb >= 100 ~ "Ultra-hypermutant")
  ) %>%
  # join with tp53_compare
  inner_join(
    select(tp53_compare,
           Kids_First_Biospecimen_ID_DNA,
           Kids_First_Biospecimen_ID_RNA,
           tp53_score),
    by = "Kids_First_Biospecimen_ID_DNA"
  ) %>%
  # rename RNA column for joining
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA) %>%
  # join with extend_scores using RNA column
  inner_join(
    select(extend_scores,
           Kids_First_Biospecimen_ID = SampleID, # rename for joining
           telo_score = NormEXTENDScores
    )
  )


# Prepare combined data for visualization
plot_df <- tp53_telo_mutator_df %>%
  # add in histology information
  inner_join(
    select(histologies_df,
           Kids_First_Biospecimen_ID,
           cancer_group,
           broad_histology)
  ) %>%
  # add in palette information
  inner_join(by = c("broad_histology", "cancer_group"),
    select(histologies_palette_df,
           cancer_group,
           cancer_group_display,
           cancer_group_hex,
           broad_histology)
  ) %>%
  # filter to cancer groups of interest
  filter(cancer_group_display %in% cancer_groups_to_plot) %>%
  # duplicate the tp53 scores column so we can eventually order can groups by it
  mutate(tp53_forordering = tp53_score) %>%
  # we want a single column for all scores so we can facet by scores
  # first, change naming for facet labels:
  rename(`Telomerase score` = telo_score,
         `TP53 score` = tp53_score) %>%
  gather(contains("score"),
         key = "score_type",
         value = "score") %>%
  # order TP53 on top
  mutate(score_type = fct_relevel(score_type, "TP53 score")) %>%
  # wrap cancer group label
  mutate(cancer_group_display = str_wrap(cancer_group_display, cg_wrap))


# Define colors to use
legend_colors <- c(Normal      = "grey40",
                   Hypermutant = "orange",
                   `Ultra-hypermutant` = "red")
# Other plot parameters which need to be re-introduced in legend:
normal_alpha <- 0.7
normal_size <- 0.75
mutator_size <- 1
point_stroke <- 0.25
normal_pch <- 19
mutator_pch <- 21
jitter_width <- 0.15 # not in legend but often in main plot

# Boxplot with overlayed jitter with colors "mapped" to `mutator`
set.seed(12) # reproducible jitter to ensure we can see all N=6 points in BOTH facets
tp53_telo_tmb_boxplot <- ggplot(plot_df) +
  aes(
    x = fct_reorder(cancer_group_display, tp53_forordering), # order cancer groups by tp53 score
    y = score
  ) +
  geom_boxplot(
    outlier.shape = NA, # no outliers
    color = "grey20",    # dark grey color
    size = 0.2
  ) +
  # Separate out jitters so that the mutant layers are ON TOP OF normal
  geom_jitter(
    data = plot_df[plot_df$mutator == "Normal",],
    width = jitter_width,
    alpha = normal_alpha,
    size = normal_size,
    pch = normal_pch,
    color = legend_colors["Normal"],
    stroke = point_stroke
  ) +
  geom_jitter(
    data = plot_df[plot_df$mutator == "Hypermutant",],
    width = jitter_width,
    pch = mutator_pch,
    size = mutator_size,
    fill = legend_colors["Hypermutant"],
    stroke = point_stroke
  ) +
  geom_jitter(
    data = plot_df[plot_df$mutator == "Ultra-hypermutant",],
    width = jitter_width,
    pch = mutator_pch,
    size = mutator_size,
    fill = legend_colors["Ultra-hypermutant"],
    stroke = point_stroke
  ) +
  labs(x = "Cancer group",
       y = "Score") +
  facet_wrap(~score_type, nrow = 2) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = rel(0.6)),
    # Sizing for compilation - small figure export
    axis.text.y = element_text(size = rel(0.6)),
    axis.title = element_text(size = rel(0.65)),
    strip.text = element_text(size = rel(0.75)),
    axis.line = element_line(size = rel(0.5)),
    axis.ticks = element_line(size = rel(0.5))
  )

# Export plot
ggsave(tp53_telomerase_scores_boxplot_pdf,
       tp53_telo_tmb_boxplot,
       width = 4.875, height = 3.25,
       # compilation:
       useDingbats = FALSE)



# Make a legend for the grey/orange/red since this was not done with normal mapping
# We have to make a "fake" plot for this to extract the legend from

# Count numbers of each category to show in legend
tp53_plot_legend_df <- plot_df %>%
  # keep 1 category only for counting
  filter(score_type == "TP53 score") %>%
  # column to use for ordering mutator levels
  mutate(mutator_order = case_when(
    mutator == "Normal" ~ 1,
    mutator == "Hypermutant" ~ 2,
    mutator == "Ultra-hypermutant" ~ 3
  )) %>%
  # count
  group_by(mutator) %>%
  mutate(mutator_count = n()) %>%
  ungroup() %>%
  # update labeling with count information
  mutate(mutator = glue::glue("{mutator} (N = {mutator_count})"),
         mutator_factor = factor(mutator),
         mutator_factor = fct_reorder(mutator_factor, mutator_order))


legend_name <- "Mutation status"
tp53_plot_for_legend <- ggplot(tp53_plot_legend_df) +
  aes(x = cancer_group, y = tmb, shape = mutator_factor, fill = mutator_factor, color = mutator_factor, size = mutator_factor, alpha = mutator_factor) +
  geom_point(size =3) +
  scale_size_manual(name = legend_name, values = c(normal_size, mutator_size, mutator_size/10))+
  scale_shape_manual(name = legend_name, values = c(normal_pch, mutator_pch, mutator_pch)) +
  scale_alpha_manual(name = legend_name, values = c(normal_alpha, 1, 1))+
  scale_color_manual(name = legend_name,values = c(unname(legend_colors["Normal"]), "black", "black")) +
  scale_fill_manual(name = legend_name, values = c("black", unname(legend_colors["Hypermutant"]), unname(legend_colors["Ultra-hypermutant"]))) +
  # theme to remove gray background. this strategy works
  theme_classic() +
  theme(
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.key.size = unit(15, "points"), 
    legend.position = "bottom"
  ) 


legend <- cowplot::get_legend(tp53_plot_for_legend)

# Export legend
pdf(tp53_telomerase_scores_boxplot_legend_pdf, width = 4, height = 0.5, useDingbats = FALSE)
print(cowplot::ggdraw(legend)) # add print() to ensure formatting
dev.off()


## Survival analysis forest plot -----------------------------------------------------------------
resect_ref <- "Tumor resection: Biopsy (ref)"
lgg_ref <- "LGG group: non-LGG (ref)"
hgg_ref <- "HGG group: non-HGG (ref)"

# Set up ordering and labels for y-axis
term_order <- rev(c("tp53_score",
                    "telomerase_score",
                    "extent_of_tumor_resectionGross/Near total resection",
                    "extent_of_tumor_resectionPartial resection",
                    resect_ref,
                    "lgg_groupLGG",
                    lgg_ref,
                    "hgg_groupHGG",
                    hgg_ref))

term_labels <- rev(c("TP53 score",
                     "Telomerase score",
                     "Tumor resection: Total",
                     "Tumor resection: Partial",
                     resect_ref,
                     "LGG group: LGG",
                     lgg_ref,
                     "HGG group: HGG",
                     hgg_ref))


# Get n and event info from glance outpout
survival_n <- broom::glance(survival_result$model) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result$model) %>%
  # add references
  add_row(term = resect_ref, estimate = 0) %>%
  add_row(term = lgg_ref, estimate = 0) %>%
  add_row(term = hgg_ref, estimate = 0) %>%
  #remove unknown resection from plot
  filter(term != "extent_of_tumor_resectionUnavailable") %>%
  mutate(estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         # significance indicator column for filling points.
         significant = case_when(p.value <= 0.05 ~ "TRUE",
                                 p.value > 0.05 ~ "FALSE",
                                 is.na(p.value) ~ "REF"),
         # y-axis factor re-labeling
         term = factor(term,
                       levels = term_order,
                       labels = term_labels)
  )



# Forest plot of the model
forest_plot <- ggplot(survival_df) +
  aes(x = estimate,
      y = term,
      fill = significant
  ) +
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.25
  ) +
  geom_point(
    size = 4.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white",
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) +
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1,
    linetype = 3, 
    size = 0.25
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) +
  # log-scale the x-axis
  scale_x_log10(
    # axis label was being cut off - this fixes it
    limits=c(5e-3, 1e2)
  ) +
  ggpubr::theme_pubr() +
  theme(
    plot.subtitle = element_text(face = "bold"), 
    # thinner axes, ticks for compilation
    axis.line = element_line(size = rel(0.25)),
    axis.ticks = element_line(size = rel(0.25))
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels.
    p_string = if_else(
      p.value >= 0.001,
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  # remove CI for refs
  mutate(value = ifelse(grepl("ref", term), NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) +
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                      ",
                      "HR (95% CI)                   P-value")
  ) +
  ggpubr::theme_pubr() +
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get covered
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  )

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5))


# Export plot
ggsave(survival_plot_pdf, forest_panels, width = 11, height = 3.5)

## Export CSVs for Zenodo upload ------------------------------

# Panel 4A: ROC curve
# no sample information so no arranging is needed
readr::write_csv(roc_df, fig4a_csv)

# Panel 4B: TP53 scores violin plots
tp53_scores_data %>%
  # reorder columns so sample id first, and rename it to match what the metadata will have
  dplyr::select(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n and N=## from tp53_altered so that CSV is properly formatted
  dplyr::mutate(tp53_altered = stringr::str_replace(tp53_altered, "\n.+", "")) %>%
  # export
  readr::write_csv(fig4b_csv)


# Panel 4C: TP53 expression plots
tp53_expression_data %>%
  # reorder columns so sample id first, and rename it to match what the metadata will have
  dplyr::select(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n and N=## from tp53_altered so that CSV is properly formatted
  dplyr::mutate(tp53_altered = stringr::str_replace(tp53_altered, "\n.+", "")) %>%
  # export
  readr::write_csv(fig4c_csv)


# Panel 4D: Box/jitter plots of TP53 and EXTEND scores across cancer groups
plot_df %>%
  # reorder columns so the RNA ID (this df also has `_DNA`) is first
  dplyr::select(Kids_First_Biospecimen_ID, everything()) %>%
  # arrange on RNA sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n from cancer_group_display so that CSV is properly formatted
  # Note that column has only wrapped the cancer group name; there is no `(N=##)` to remove
  dplyr::mutate(cancer_group_display = stringr::str_replace(cancer_group_display, "\n", " ")) %>%
  # export
  readr::write_csv(fig4d_csv)


# Panel 4F: TP53 and telomerase forest plot
survival_result$data %>%
  # order columns and only keep columns that are in the model: 
  # OS_status and OS_years
  # tp53_score+telomerase_score+extent_of_tumor_resection+lgg_group+hgg_group"
  dplyr::select(Kids_First_Biospecimen_ID, 
                OS_years, 
                OS_status,
                tp53_score, 
                telomerase_score, 
                extent_of_tumor_resection,
                lgg_group,
                hgg_group) %>%
  # recode OS_status back to LIVING/DECEASED
  dplyr::mutate(OS_status = ifelse(
    OS_status == 1, 
    "LIVING", 
    "DECEASED")
  ) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig4f_csv)
