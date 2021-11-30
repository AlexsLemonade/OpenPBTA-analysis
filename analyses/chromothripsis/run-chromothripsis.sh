#!/bin/bash

# Module author: Laura Egolf
# Bash script author: Jaclyn Taroni

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Process Manta output
Rscript --vanilla 01-process-sv-file.R

# Shatterseek
Rscript --vanilla 02-run-shatterseek-and-classify-confidence.R

# Plotting - by disease label
Rscript --vanilla -e "rmarkdown::render('03-plot-chromothripsis-by-histology.Rmd', clean = TRUE)"

# Plotting - breakpoint data
Rscript --vanilla -e "rmarkdown::render('04-plot-chromothripsis-and-breakpoint-data.Rmd', clean = TRUE)"
