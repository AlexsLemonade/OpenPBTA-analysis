#!/bin/bash

# Run Jin

# Generate TSNE figures 

set -e
set -o pipefail

# Set the working directory to the directory of this file
 cd "$(dirname "${BASH_SOURCE[0]}")"

# Plot TSNE for cancer group of interest 
Rscript -e "rmarkdown::render('tsne_figures.Rmd')"
