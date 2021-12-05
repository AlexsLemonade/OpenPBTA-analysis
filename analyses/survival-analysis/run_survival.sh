#!/bin/bash

# K S Gaonkar, J Rokita

# Run survival-analysis

set -e
set -o pipefail

# Plot survial curves of subtypes in HGG/DMG samples
Rscript -e "rmarkdown::render('analyses/survival-analysis/survival-analysis_HGG_DMG.Rmd')"
Rscript -e "rmarkdown::render('analyses/survival-analysis/survival-analysis_histologies.Rmd')"

