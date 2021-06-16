#!/bin/bash

# J. Taroni for ALSF CCDL 
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# In CI we'll run an abbreviated version of the de novo signatures extraction
ABBREVIATED_MUTSIGS=${OPENPBTA_QUICK_MUTSIGS:-0}


# Run the mutational signatures analysis using existing signatures
Rscript -e "rmarkdown::render('01-known_signatures.Rmd', clean = TRUE)"

# Split up the consensus MAF files by experimental strategy (writes to scratch)
Rscript --vanilla 02-split_experimental_strategy.R

# Run the shell script for determining the number of signatures to use
# with a low number of iterations if run in CI
# argument 0 --> run denovo goodness of fit
QUICK_MUTSIGS=$ABBREVIATED_MUTSIGS ANALYSIS=0 bash 03-de_novo_range_of_nsignatures.sh

# Run the shell script to perform de novo extraction 
# with a low number of iterations if run in CI
# argument 1 --> run more robust denovo extraction
QUICK_MUTSIGS=$ABBREVIATED_MUTSIGS ANALYSIS=1 bash 03-de_novo_range_of_nsignatures.sh

# Process results from de novo extraction 
Rscript -e "rmarkdown::render('04-analyze_de_novo.Rmd', clean = TRUE, params = list(is_ci = ${ABBREVIATED_MUTSIGS}))"


# Next steps: Fitting the 8 known CNS signatures