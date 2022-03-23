# Stephanie J. Spielman for CCDL, 2022
# This script creates panels for Figure S2
# https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1254

#A: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/pbta-caller-comparison/pbta-vaf_cor_matrix.png
#B: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/pbta-caller-comparison/pbta-vaf_distribution_plot.png
#C: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/pbta-caller-comparison/pbta-upset-plot.png
#D: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/tcga-caller-comparison/tcga-vaf-cor-matrix.png
#E: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/tcga-caller-comparison/tcga-vaf-distribution-plot.png
#F: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/plots/tcga-caller-comparison/tcga-upset-plot.png
#G&H: (split into two panels!) https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/tmb-compare/plots/tmb-cdf-pbta-tcga.png


# Load libraries ---------------------
library(tidyverse)


# Set up paths and file names -----------------------

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# make output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "supp_figure2_panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses")

snv_callers_dir <- file.path(analyses_dir, "snv-callers")


# Load files and source relevant utils ------------------------




