#!/bin/bash
# JN Taroni for ALSF CCDL 2022
# Create tables of git contributions to the current branch

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Shell script wrapped around git shortlog and git log
bash 01-count-contributions.sh

# Run notebook with data wrangling
Rscript -e "rmarkdown::render('02-format-contributions.Rmd', clean = TRUE)"

# Run notebook reordering author list to account for git contributions
Rscript -e "rmarkdown::render('03-set-authorship-order.Rmd', clean = TRUE)"

# Run notebook generating author list as TSV for submission
Rscript -e "rmarkdown::render('04-get-author-information.Rmd', clean = TRUE)"

