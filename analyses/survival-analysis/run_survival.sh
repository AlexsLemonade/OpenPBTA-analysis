#!/bin/bash

# Jo Lynne Rokita and Run Jin

# Run survival analysis 

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Plot survial curves of subtypes in HGG/DMG samples
Rscript -e "rmarkdown::render('survival-analysis_HGG_DMG.Rmd')"

# Plot survival curves of histology
Rscript -e "rmarkdown::render('survival-analysis_histology.Rmd')"

# Plot survival curves using immune results 
Rscript -e "rmarkdown::render('survival-analysis_immune.Rmd')"
