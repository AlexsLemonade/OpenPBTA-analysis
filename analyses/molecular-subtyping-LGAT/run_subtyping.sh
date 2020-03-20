#!/bin/bash

set -e
set -o pipefail

# Run notebook to get molecular subtype for LGAT samples 

Rscript -e "rmarkdown::render('analyses/molecular-subtyping-LGAT/01-make-lgat-final-table.Rmd')"

