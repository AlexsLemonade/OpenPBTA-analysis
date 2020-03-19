#!/bin/bash

set -e
set -o pipefail

# Run notebook to get reclassified metadata 

Rscript -e "rmarkdown::render('analyses/molecular_classifying_EWS/01-reclassify_as_ewings.Rmd')"


