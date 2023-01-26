#!/bin/bash

# SJ Spielman (CCDL) 2022-2023

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript -e "rmarkdown::render('01_explore-tumor-purity.Rmd', clean = TRUE)"

Rscript -e "rmarkdown::render('02_explore-tumor-purity-thresholds.Rmd', clean = TRUE)"
