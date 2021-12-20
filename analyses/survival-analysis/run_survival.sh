#!/bin/bash

# K S Gaonkar

# Run fusion_filtering

set -e
set -o pipefail

# Plot survial curves of subtypes in HGG/DMG samples
Rscript -e "rmarkdown::render('survival-analysis_HGG_DMG.Rmd')"

Rscript -e "rmarkdown::render('survival-analysis_histology.Rmd')"
