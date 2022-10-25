#!/bin/bash

# Bethell and Taroni for CCDL 2019
# Run the dimension reduction plotting pipeline specifically in CI

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run stranded RSEM file with low perplexity
Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead rsem_stranded \
  --output_directory results \
  --perplexity 3

# Run poly-A kallisto file, skipping t-SNE
Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.polya.rds \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead kallisto_polyA \
  --output_directory results \
  --skip_tsne \
  --neighbors 2

# generate plot lists for both cases above
Rscript --vanilla scripts/get-plot-list.R  \
  --input_directory results \
  --filename_lead rsem_stranded \
  --output_directory plots/plot_data \
  --color_variable broad_histology

Rscript --vanilla scripts/get-plot-list.R  \
  --input_directory results \
  --filename_lead kallisto_polyA \
  --output_directory plots/plot_data \
  --color_variable broad_histology

# Run all the multipanel plots!
bash 03-multipanel-plots.sh

# Exploration of batch effects
Rscript --vanilla -e 'rmarkdown::render("04-explore-sequencing-center-effects.Rmd")'
