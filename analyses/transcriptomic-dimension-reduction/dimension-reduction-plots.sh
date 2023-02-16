#!/bin/bash

# Bethell and Taroni for CCDL 2019
# Run the dimension reduction plotting pipeline. Samples will be colored by
# user-specified variable. It will be broad_histology by default.

set -e
set -o pipefail

COLORVAR=${COLOR:-broad_histology}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# PCA, UMAP, t-SNE step
bash 01-dimension-reduction.sh

# Generate plot lists to be used to make multipanel plots
COLOR=${COLORVAR} bash 02-get-dimension-reduction-plot-lists.sh
# Make multipanel plots and save as PDFs
bash 03-multipanel-plots.sh

# Exploration of batch effects
Rscript --vanilla -e 'rmarkdown::render("04-explore-sequencing-center-effects.Rmd")'

# Exploration of UMAPs if mitochondrial genes are removed
Rscript --vanilla -e 'rmarkdown::render("05-seq-center-mitochondrial-genes.Rmd")'

# UMAP visualization with the tumor purity filtered dataset
Rscript --vanilla -e 'rmarkdown::render("06-umap-tumor-purity-threshold.Rmd")'

