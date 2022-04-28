#!/bin/bash

# Jo Lynne Rokita (D3b), Run Jin (D3b), and Stephanie Spielman (CCDL)

# Run survival analysis 

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Plot survial curves of subtypes in HGG/DMG samples
Rscript -e "rmarkdown::render('survival-analysis_subtypes.Rmd', clean = TRUE)"

# Build survival models to assess tp53, telomerase, cancer group, and HGG group effects
Rscript -e "rmarkdown::render('survival-analysis_tp53_telomerase.Rmd', clean = TRUE)"

# Build survival models to assess molecular subtype, PDL1 expression, and immune cell fraction effects
Rscript -e "rmarkdown::render('survival-analysis_immune.Rmd', clean = TRUE)"
