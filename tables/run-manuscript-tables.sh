#!/bin/bash
# Module authors: Run Jin (D3b), Stephanie J. Spielman (ALSF CCDL), and Jo Lynne Rokita (D3b)
# Shell script author: Jo Lynne Rokita (D3b)
# 2022

# This script runs the steps for generating manuscript tables.

set -e
set -o pipefail

# run the notebook to investigate hypermutator BS_F0GNWEJJ
Rscript -e "rmarkdown::render('util/BS_F0GNWEJJ_genomic_investigation.Rmd')"

# run the notebook to create manuscript tables
Rscript -e "rmarkdown::render('output_tables.Rmd')"
