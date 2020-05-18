#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run notebook to get reclassified metadata 

Rscript -e "rmarkdown::render('01-reclassify_as_ewings.Rmd')"


