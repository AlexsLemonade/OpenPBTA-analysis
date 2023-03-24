# S. Spielman for ALSF CCDL 2023
#
# Makes pdf panels for TP53 and telomerase scores across cancer groups, focusing only
#  on high tumor purity samples

library(tidyverse)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs7", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directories
analyses_dir <- file.path(root_dir, "analyses")
tp53_dir <- file.path(analyses_dir, "tp53_nf1_score")
telomerase_dir <- file.path(analyses_dir, "telomerase-activity-prediction")

# Palette directory
palette_dir <- file.path(root_dir, "figures", "palettes")

# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS7g_csv <- file.path(zenodo_tables_dir, "figure-S7g-data.csv")

# Output PDF filenames
output_pdf <- file.path(output_dir, "tp53_telomerase_boxplots_tumor-purity-threshold.pdf")
output_legend_pdf <- file.path(output_dir, "tp53_telomerase_boxplots_legend_tumor-purity-threshold.pdf")


# Read in clinical data and associated palette
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)
histologies_palette_df <- read_tsv(file.path(palette_dir, "broad_histology_cancer_group_palette.tsv"))


# Read in tp53 data file
tp53_df <- read_tsv(file.path(tp53_dir,
                              "results",
                              "tumor-purity-threshold",
                              "tp53_altered_status_tumor-purity-threshold.tsv"))

# Read in EXTEND scores
extend_df <- read_tsv(file.path(telomerase_dir,
                                "results",
                                "TelomeraseScores_PTBAStranded_FPKM_thresholded.txt"))

# Read in TMB for highlighting points in boxplots
tmb_coding_df <- read_tsv(file.path(data_dir, "pbta-snv-consensus-mutation-tmb-coding.tsv"))

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
  # note that there should only be 2 Hypermutators after the filtering on tumor purity
  mutate(
    mutator = case_when(
      tmb < 10 ~ "Normal",
      tmb >= 10 & tmb < 100 ~ "Hypermutant",
      tmb >= 100 ~ "Ultra-hypermutant")
  ) %>%
  # join with tp53_compare
  inner_join(
    select(tp53_df,
           Kids_First_Biospecimen_ID_DNA,
           Kids_First_Biospecimen_ID_RNA,
           tp53_score),
    by = "Kids_First_Biospecimen_ID_DNA"
  ) %>%
  # rename RNA column for joining
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA) %>%
  # join with extend_scores using RNA column
  inner_join(
    select(extend_df,
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
           cancer_group)
  ) %>%
  # add in palette information
  inner_join(
    select(histologies_palette_df,
           cancer_group,
           cancer_group_display)
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
                   Hypermutant = "orange")
# Other plot parameters which need to be re-introduced in legend:
normal_alpha <- 0.7
normal_size <- 0.75
mutator_size <- 1
point_stroke <- 0.25
normal_pch <- 19
mutator_pch <- 21
jitter_width <- 0.15 # not in legend but often in main plot

# Boxplot with overlayed jitter with colors "mapped" to `mutator`
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
  labs(x = "Cancer group",
       y = "Score") +
  facet_wrap(~score_type, nrow = 2) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.5)),
    # Sizing for compilation - small figure export
    axis.text.y = element_text(size = rel(0.5)),
    axis.title = element_text(size = rel(0.5)),
    strip.text = element_text(size = rel(0.5)),
    axis.line = element_line(size = rel(0.5)),
    axis.ticks = element_line(size = rel(0.5))
  )

# Export plot
ggsave(output_pdf,
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
    mutator == "Hypermutant" ~ 2
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
  geom_point(size = 3) +
  scale_size_manual(name = legend_name, values = c(normal_size, mutator_size)) +
  scale_shape_manual(name = legend_name, values = c(normal_pch, mutator_pch)) +
  scale_alpha_manual(name = legend_name, values = c(normal_alpha, 1)) +
  scale_color_manual(name = legend_name,values = c(unname(legend_colors["Normal"]), "black")) +
  scale_fill_manual(name = legend_name, values = c("black", unname(legend_colors["Hypermutant"]))) +
  # theme to remove gray background. this strategy works
  theme_classic() +
  theme(
  # Sizing for compilation - small figure export
    legend.text = element_text(size = rel(0.45)),
    legend.title = element_text(size = rel(0.6)),
    legend.key.size = unit(10, "points")
)


legend <- cowplot::get_legend(tp53_plot_for_legend)

# Export legend
pdf(output_legend_pdf, width = 1.3, height = 0.8, useDingbats = FALSE)
cowplot::ggdraw(legend)
dev.off()

# Export CSV for Zenodo upload
# no samples so nothing to arrange
plot_df %>%
  # remove \n
  dplyr::mutate(cancer_group_display = stringr::str_replace(cancer_group_display, "\n", " ")) %>%
  # arrange on RNA ID column, but bring both to front
  dplyr::select(Kids_First_Biospecimen_ID, Kids_First_Biospecimen_ID_DNA, everything()) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(figS7g_csv)


