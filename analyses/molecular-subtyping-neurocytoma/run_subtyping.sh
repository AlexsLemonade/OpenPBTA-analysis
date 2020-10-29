#!/bin/bash

set -e
set -o pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Run notebook to get molecular subtype for Neurocytoma samples 
Rscript -e "rmarkdown::render('01-neurocytoma-subtyping.Rmd')"
