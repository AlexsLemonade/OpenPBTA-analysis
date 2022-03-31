# S. Spielman for CCDL 2022
#
# Makes pdf panel of TP53 scores with (ultra)hypermutators emphasized for inclusion in Figure 4

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

# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

# PDF panels
tp53_pdf <- file.path(output_dir, "tp53_panel.pdf")
tp53_legend_pdf <- file.path(output_dir, "tp53_legend.pdf")


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



#### Read in and prepare data --------------------------------------------------
# define paths and read in TP53 scores
tp53_dir <- file.path(
  analyses_dir,
  "tp53_nf1_score"
)
tp53_scores_file <- file.path(
  tp53_dir,
  "results",
  "tp53_altered_status.tsv"
)
tp53_scores_df <- read_tsv(tp53_scores_file)


# define paths and read in TMB information for identifying (ultra)hypermutants. We define hyper as > 10 and ultrahyper as > 100 (https://doi.org/10.1016/j.cell.2017.09.048)
# we focus on coding TMB here (nothing is >10, let alone >100, in all mutations)
tmb_coding_file <- file.path(
  data_dir, 
  "pbta-snv-consensus-mutation-tmb-coding.tsv"
)
tmb_coding_df <- read_tsv(tmb_coding_file)


# Rename DNA column for joining purposes and select columns of interest
tp53_scores_df_subset <- tp53_scores_df %>%
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  select(Kids_First_Biospecimen_ID, 
         tp53_score) %>%
  distinct() %>%
  # remove all NAs, which exist in IDs and tp53
  drop_na()
  
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
      tmb >= 100 ~ "Ultrahypermutant")
  )

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
  mutate(cancer_group = str_wrap(cancer_group, 20))


#### Plot and associated legend --------------------------------------------------


# Define colors to use
legend_colors <- c(Normal      = "grey40",
                   Hypermutant = "orange",
                   Ultrahypermutant = "red")

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
ggsave(tp53_pdf,
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
pdf(tp53_legend_pdf, width = 6, height = 3)
cowplot::ggdraw(legend)
dev.off()







