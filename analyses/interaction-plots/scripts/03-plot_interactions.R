# Create plot of co-occurence and mutual exclusivity
#
# JA Shapiro for ALSF - CCDL
#
# 2019
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
# Rscript analyses/interaction-plots/02-plot_interactions.R \
#   --infile analysis/interaction-plots/results/cooccur.tsv \
#   --outfile analysis/interaction-plots/results/cooccur.png

#### Initial Set Up
# Establish base dir

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)

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
    help = "File path where geneXdisease table is located (optional)",
    metavar = "character"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

cooccur_file <- opts$infile
plot_file <- opts$outfile

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
# order genes the same way, in case we ant to use those
genes <- stringr::str_extract(labels, "^.+?\\b")
genes <- genes[order(label_counts, decreasing = TRUE)]

cooccur_df <- cooccur_df %>%
  dplyr::mutate(
    gene1 = factor(gene1, levels = genes),
    gene2 = factor(gene2, levels = genes),
    label1 = factor(label1, levels = labels),
    label2 = factor(label2, levels = labels)
  )


# create scales for consistent sizing
xscale <- cooccur_df$label1 %>%
  as.character() %>%
  unique() %>%
  c(1:(opts$plotsize - length(.)))
yscale <- cooccur_df$label2 %>%
  as.character() %>%
  unique() %>%
  c(1:(opts$plotsize - length(.)), .)

### make plot
cooccur_plot <- ggplot(
  cooccur_df,
  aes(x = label1, y = label2, fill = cooccur_score)
) +
  geom_tile(color = "white", size = 1) +
  scale_x_discrete(
    position = "top",
    limits = xscale,
    breaks = unique(cooccur_df$label1)
  ) + # blank unused sections.
  scale_y_discrete(
    limits = yscale,
    breaks = unique(cooccur_df$label2)
  ) +
  scale_fill_distiller(
    type = "div",
    palette = 5,
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

disease_file = opts$disease_table
disease_df <-
  readr::read_tsv(disease_file, col_types = readr::cols()) %>%
  dplyr::mutate(gene = factor(gene, levels = genes))

disease_plot <- ggplot(
  disease_df,
  aes(x = gene, y = mutant_samples, fill = disease)) + 
  geom_col() +
  theme_classic()


