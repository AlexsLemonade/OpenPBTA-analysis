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

# Run the shell script that is for determining the number of signatures to use
# with a low number of iterations
QUICK_MUTSIGS=$ABBREVIATED_MUTSIGS bash 03-de_novo_range_of_nsignatures.sh
