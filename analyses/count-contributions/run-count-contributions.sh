#!/bin/bash
# JN Taroni for ALSF CCDL 2022
# Create tables of git contributions to the current branch

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Shell script wrapped around git shortlog and git log
bash 01-count-contributions.sh

ls ../../scratch/count-contributions

# Notebook that cleans things up
# Rscript -e "rmarkdown::render('02-format-contributions.Rmd', clean = TRUE)"
