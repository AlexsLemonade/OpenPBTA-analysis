#!/bin/bash

# K S Gaonkar, J Rokita

# Run survival-analysis

set -e
set -o pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Plot survial curves of subtypes in HGG/DMG samples
Rscript -e "rmarkdown::render('survival-analysis_HGG_DMG.Rmd')"

# Plot survival curves based on TP53 scores and telomerase scores
Rscript -e "rmarkdown::render('survival-analysis_histologies.Rmd')"

