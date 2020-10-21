#!/bin/bash

# Josh Shapiro for CCDL 2019
#
# Runs 01-generate-independent-specimens.R with default settings.

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript -e "rmarkdown::render('00-repeated-samples.Rmd', clean = TRUE)"

Rscript 01-generate-independent-specimens.R \
  -f ../../data/pbta-histologies.tsv \
  -o results

# adding indepedent list of rnaseq samples
Rscript 02-generate-independent-rnaseq.R \
  --histology_file ../../data/pbta-histologies.tsv \
  --output_directory results \
  --independent_dna_sample_df ../../data/independent-specimens.wgswxs.primary-plus.tsv  

