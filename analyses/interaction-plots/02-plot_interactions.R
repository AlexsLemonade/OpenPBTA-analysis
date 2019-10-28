# Create plot of co-occurence and mutual exclusivity
#
# JA Shapiro for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --infile
#                 
#
# Command line example:
#
# Rscript analyses/interaction-plots/02-plot_interactions.R \
#   --infile 

#### Initial Set Up
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
script_root <- file.path(root_dir, "analyses", "interaction-plots")

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)

option_list <- list(
  make_option(
    opt_str = "--infile", type = "character", 
    default = file.path("analyses", "interaction-plots", "results", "cooccurance.tsv"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            where cooccurence summary table is located",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outfile", type = "character", 
    default = file.path("analyses", "interaction-plots", "results", "cooccurance.png"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            where output plot will be located. Extension specifies format of plot",
    metavar = "character"
  ),
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

cooccur_file <- file.path(root_dir, opts$infile)
plot_file <- file.path(root_dir, opts$infile)

cooccur_df <- readr::read_tsv(cooccur_file)
  
### make plot
ggplot(cooccurence_df, aes(x = gene1, y = gene2, fill = cooccur_score))+
  geom_tile(color = "white", size = 1) +
  xlab('') + 
  ylab('') +
  scale_x_discrete(position = "top") +
  scale_fill_distiller(type = "div", palette = 5, 
                       limits = c(-20, 20),
                       oob = scales::squish) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        legend.justification=c(1,0), 
        legend.position=c(1,0),
        legend.key.size = unit(2, "char"))

ggsave(filename = plot_file)