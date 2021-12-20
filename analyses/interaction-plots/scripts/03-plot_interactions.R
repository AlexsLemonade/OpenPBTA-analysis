# Create plot of co-occurence and mutual exclusivity
#
# JA Shapiro for ALSF - CCDL
#
# 2019-2020
#
# Option descriptions
#
# --infile The input file  with a summaries of gene-gene mutation co-occurence,
#    minimally including gene1, gene2, and cooccur_score columns.
#
# --outfile The output plot location. Specify type of file with the extension
#   (.png or .pdf, most likely).
#
# --plotsize The number of rows and columns in the expected plot, for scaling.
#   Larger numbers will create smaller boxes for the heatmap tiles.
#
# Command line example:
#
# Rscript analyses/interaction-plots/03-plot_interactions.R \
#   --infile analysis/interaction-plots/results/cooccur.tsv \
#   --outfile analysis/interaction-plots/results/cooccur.png

#### Initial Set Up

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)
library(patchwork)

# define options
option_list <- list(
  make_option(
    opt_str = "--infile",
    type = "character",
    help = "File path where cooccurence summary table is located",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outfile",
    type = "character",
    help = "File path where output plot will be located. Extension specifies format of plot",
    metavar = "character"
  ),
  make_option(
    opt_str = "--plotsize",
    default = "50",
    type = "numeric",
    help = "Relative size of plots; number of rows and columns to be plotted",
    metavar = "character"
  ),
  make_option(
    opt_str = "--disease_table",
    type = "character",
    default = NA,
    help = "File path where gene X disease table is located (optional)",
    metavar = "character"
  ),
  make_option(
    opt_str = "--disease_plot",
    type = "character",
    default = NA,
    help = "File path where gene X disease plot should be placed (required if --disease_table specified)",
    metavar = "character"
  ),
  make_option(
    opt_str = "--combined_plot",
    type = "character",
    default = NA,
    help = "File path where gene X disease plot should be placed (required if --disease_table specified)",
    metavar = "character"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opts$disease_table)){
  if (is.na(opts$disease_plot) | is.na(opts$combined_plot)){
    stop("If disease_table is specified, disease_plot and/or combined plot must also be specified")
  }
}

cooccur_file <- opts$infile
plot_file <- opts$outfile

# get root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


cooccur_df <-
  readr::read_tsv(cooccur_file, col_types = readr::cols()) %>%
  dplyr::mutate(
    mut1 = mut11 + mut10,
    mut2 = mut11 + mut01,
    label1 = paste0(gene1, " (", mut1, ")"),
    label2 = paste0(gene2, " (", mut2, ")")
  )

labels <- unique(c(cooccur_df$label1, cooccur_df$label2))

# check the order of the labels to be decreasing by mut count
label_counts <- as.numeric(stringr::str_extract(labels, "\\b\\d+\\b"))
labels <- labels[order(label_counts, decreasing = TRUE)]
# order genes the same way, in case we want to use those
genes <- stringr::str_extract(labels, "^.+?\\b")
genes <- genes[order(label_counts, decreasing = TRUE)]

cooccur_df <- cooccur_df %>%
  dplyr::mutate(
    gene1 = factor(gene1, levels = genes),
    gene2 = factor(gene2, levels = genes),
    label1 = factor(label1, levels = labels),
    label2 = factor(label2, levels = labels)
  )


# Get color palettes

palette_dir <- file.path(root_dir, "figures", "palettes")
divergent_palette <- readr::read_tsv(file.path(palette_dir, "divergent_color_palette.tsv"),
                                     col_types = readr::cols())
divergent_colors <- divergent_palette %>%
  dplyr::filter(color_names != "na_color") %>%
  dplyr::pull(hex_codes)
na_color <- divergent_palette %>%
  dplyr::filter(color_names == "na_color") %>%
  dplyr::pull(hex_codes)

histologies_color_key_df <- readr::read_tsv(file.path(palette_dir,
                                                      "broad_histology_cancer_group_palette.tsv"),
                                            col_types = readr::cols())

# create scales for consistent sizing
# The scales will need to have opts$plotsize elements,
# so after getting the unique list, we concatenate on extra elements.
# for convenience, these are just numbers 1:n
# where n is the number of extra labels needed for the scale
xscale <- cooccur_df$label1 %>%
  as.character() %>%
  unique() %>%
  c(1:(opts$plotsize - length(.)))
yscale <- cooccur_df$label2 %>%
  as.character() %>%
  unique() %>%
  # the concatenated labels need to be at the front of the Y scale,
  # since this will be at the bottom in the plot.
  c(1:(opts$plotsize - length(.)), .)

### make plot
cooccur_plot <- ggplot(
  cooccur_df,
  aes(x = label1, y = label2, fill = cooccur_score)
) +
  geom_tile(width = 0.7, height = 0.7) +
  scale_x_discrete(
    position = "top",
    limits = xscale,
    breaks = unique(cooccur_df$label1)
  ) + # blank unused sections.
  scale_y_discrete(
    limits = yscale,
    breaks = unique(cooccur_df$label2)
  ) +
  scale_fill_gradientn(
    colors = divergent_colors,
    na.value = na_color,
    limits = c(-10, 10),
    oob = scales::squish,
  ) +
  labs(
    x = "",
    y = "",
    fill = "Co-occurence\nscore"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(
      angle = -90,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 6),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0),
    legend.key.size = unit(2, "char")
  )

ggsave(cooccur_plot, filename = plot_file)

# if we don't have a disease table, quit
if (is.na(opts$disease_table)) {
 quit()
}
# otherwise make a gene by disease stacked bar chart

disease_file <- opts$disease_table
disease_df <-
  readr::read_tsv(disease_file, col_types = readr::cols()) %>%
  dplyr::mutate(gene = factor(gene, levels = genes))

# Previously calculated top 10
# display_diseases <- disease_df %>%
#   dplyr::select(disease, mutant_samples) %>%
#   dplyr::arrange(desc(mutant_samples)) %>%
#   dplyr::select(disease) %>%
#   unique() %>%
#   head(10) %>% # top 10 diseases with highest mutated samples
#   dplyr::pull(disease)

# Now a manual version
display_diseases <- c(
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

disease_df <- disease_df %>%
  # remove disease == NA, these are samples where
  # harmonized_diagnosis is Benign tumor, Dysplasia/Gliosis
  dplyr::filter(!is.na(disease)) %>%
  dplyr::mutate(disease_factor =
           forcats::fct_other(disease, keep = display_diseases) %>%
           forcats::fct_relevel(display_diseases)
  ) %>%
  # If you are to outline the stacked bars in anyway, all Other samples need to
  # be summarized
  dplyr::group_by(gene, disease_factor) %>%
  dplyr::summarize(mutant_samples = sum(mutant_samples)) %>%
  dplyr::ungroup()

histologies_color_key <- histologies_color_key_df$cancer_group_hex
names(histologies_color_key) <- histologies_color_key_df$cancer_group
# Adding gray for Other histologies outside the top 10 above
histologies_color_key <- c(histologies_color_key, "Other" = "#d3d3d3")


# get scale to match cooccurence plot
# Extra scale units for the case where there are fewer genes than opts$plotsize
xscale2 <- levels(disease_df$gene) %>%
  c(rep("", opts$plotsize - length(.)))

disease_plot <- ggplot(
  disease_df,
  aes(x = gene,
      y = mutant_samples,
      fill = disease_factor)) +
  geom_col(width = 0.7,
           color = "#666666") +
  labs(
    x = "",
    y = "Samples with mutations",
    fill = "Cancer Group"
  ) +
  scale_fill_manual(values = histologies_color_key) +
  scale_x_discrete(
    limits = xscale2,
    breaks = disease_df$gene
  ) +
  scale_y_continuous(expand = c(0, 0.5, 0.1, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    legend.key.size = unit(1.5, "char"))

if (!is.na(opts$disease_plot)){
  ggsave(opts$disease_plot, disease_plot)
}

# only proceed if we want a combined plot
if (is.na(opts$combined_plot)){
  quit()
}


# Modify cooccur plot to drop counts and X axis

# labels for y axis will be gene names, with extra spaces (at bottom) blank
ylabels  <- cooccur_df$gene2%>%
  as.character() %>%
  unique() %>%
  c(rep("", opts$plotsize - length(.)), .)

cooccur_plot2 <- cooccur_plot +
  scale_x_discrete(
    limits = xscale,
    breaks = c()
  ) +
  scale_y_discrete(
    limits = yscale,
    labels = ylabels
  ) +
  theme(
    plot.margin = unit(c(-3.5,0,0,0), "char") # negative top margin to move plots together
  )

# Move labels and themes for disease plot
disease_plot2 <- disease_plot +
  theme(
    axis.text.x = element_text(
      angle = -90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title.y = element_text(
      vjust = -10 # keep the label close when combined
    )
  )

# Combine plots with <patchwork>
# Layout of the two plots will be one over the other (1 column),
# with the upper plot 3/4 the height of the lower plot
combined_plot <- disease_plot2 + cooccur_plot2 +
  plot_layout(ncol = 1, heights = c(3, 4)) +
  theme( # add uniform labels
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )


ggsave(combined_plot,
       filename = opts$combined_plot,
       width = 8,
       height = 14)
