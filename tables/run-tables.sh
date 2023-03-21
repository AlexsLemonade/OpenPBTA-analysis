#!/bin/bash
# Module authors: Run Jin (D3b), Stephanie J. Spielman (ALSF CCDL), and Jo Lynne Rokita (D3b)
# Shell script author: Jo Lynne Rokita (D3b)
# 2022-2023

# This script runs the steps for generating manuscript tables.

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# testing parameter
OPENPBTA_TESTING=${OPENPBTA_TESTING:-0}

# run the notebook to investigate hypermutator BS_F0GNWEJJ
Rscript -e "rmarkdown::render(file.path('util', 'BS_F0GNWEJJ_genomic_investigation.Rmd'))"

# run the notebook to create manuscript tables, with param if testing
if [ ${OPENPBTA_TESTING} -eq 1 ]; then
    Rscript -e "rmarkdown::render('write-manuscript-tables.Rmd', params = list(release = 'testing'), clean = TRUE)"
else
    Rscript -e "rmarkdown::render('write-manuscript-tables.Rmd', clean = TRUE)"
    # Forthcoming: Copy zenodo tables if NOT in CI
fi


