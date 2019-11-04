# Create plot of co-occurence and mutual exclusivity
#
# JA Shapiro for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --infile The input file location (relative to root) with a summaries of
#   gene-gene mutation co-occurence, minimally including gene1, gene2, and
#   cooccur_score columns.
#
# --outfile The output plot location. Specify type of file with the extension
#   (.png or .pdf, most likely).
#
# Command line example:
#
# Rscript analyses/interaction-plots/02-plot_interactions.R \
#   --infile analysis/interaction-plots/results/cooccur.tsv \
#   --outfile analysis/interaction-plots/results/cooccur.png

#### Initial Set Up
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)

option_list <- list(
  make_option(
    opt_str = "--infile",
    type = "character",
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            where cooccurence summary table is located",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outfile",
    type = "character",
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            where output plot will be located. Extension specifies format of plot",
    metavar = "character"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

cooccur_file <- file.path(root_dir, opts$infile)
plot_file <- file.path(root_dir, opts$outfile)

cooccur_df <-
  readr::read_tsv(cooccur_file, col_types = readr::cols()) %>%
  dplyr::mutate(
         mut1 = mut11 + mut10,
         mut2 = mut11 + mut01,
         label1 = paste0(gene1, " (", mut1, ")"),
         label2 = paste0(gene2, " (", mut2, ")")
  )

labels <- unique(c(cooccur_df$label1, cooccur_df$label2))

cooccur_df <- cooccur_df %>%
  dplyr::mutate(
    label1 = factor(label1, levels = labels),
    label2 = factor(label2, levels = labels)
  )


### make plot
cooccur_plot <- ggplot(
  cooccur_df,
  aes(x = label1, y = label2, fill = cooccur_score)
) +
  geom_tile(color = "white", size = 1) +

  scale_x_discrete(position = "top") +
  scale_fill_distiller(
    type = "div",
    palette = 5,
    limits = c(-20, 20),
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
