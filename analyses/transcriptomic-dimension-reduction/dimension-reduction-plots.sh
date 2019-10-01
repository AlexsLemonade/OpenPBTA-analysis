#!/bin/bash

# Bethell and Taroni for CCDL 2019
# Run the dimension reduction plotting pipeline. Samples will be colored by
# user-specified variable. It will be broad_histology by default.

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
