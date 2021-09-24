#!/bin/bash

# Josh Shapiro for CCDL 2019
#
# Runs 01-generate-independent-specimens.R with default settings.
# Takes one environment variable, `OPENPBTA_BASE_SUBTYPING`, if value is 1 then
# uses histologies-base.tsv for subtyping if value is 0 runs all modules with histologies.tsv(Default)

set -e
set -o pipefail

RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# run initial script
Rscript -e "rmarkdown::render('00-repeated-samples.Rmd',params=list(base_run = ${RUN_FOR_SUBTYPING}), clean = TRUE)"

# run script to generate wgs-only lists
Rscript 01-generate-independent-specimens-wgs-only.R

# run script to generate wgs-preferred lists
Rscript 01-generate-independent-specimens-wgs-preferred.R

# run script to generate wxs-preferred lists
Rscript 01-generate-independent-specimens-wxs-preferred.R 

# run script to generate rnaseq lists
Rscript 02-generate-independent-rnaseq.R
  
# run summary on output files
Rscript -e "rmarkdown::render('03-qc-independent-samples.Rmd', clean = TRUE)"

