#!/bin/bash

set -e
set -o pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# This option controls whether on not the step that generates the LGAT only mutation 
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

if [ "$SUBSET" -gt "0" ]; then
  Rscript --vanilla 00-subset-files-for-LGAT.R
fi

# Run notebook to get molecular subtype for LGAT samples 
Rscript -e "rmarkdown::render('01-make-lgat-final-table.Rmd')"
