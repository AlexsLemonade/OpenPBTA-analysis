#!/bin/bash

# Josh Shapiro for CCDL 2019
#
# Runs 01-generate-independent-specimens.R with default settings.

set -e
set -o pipefail

Rscript -e "rmarkdown::render('analyses/independent-samples/00-repeated-samples.Rmd', clean = TRUE)"

Rscript analyses/independent-samples/01-generate-independent-specimens.R \
  -f data/pbta-histologies.tsv \
  -o analyses/independent-samples/results