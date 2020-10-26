#!/bin/bash

set -e
set -o pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Run notebook
Rscript -e "rmarkdown::render('00-craniopharyngiomas-molecular-subtype.Rmd')"
