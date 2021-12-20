#!/bin/bash

# J. Taroni for ALSF CCDL 
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Only fit for CNS signatures, which we'll do by default
# This does require us to split the MAF files by strategy, but that's all
CNS_FIT_ONLY=${OPENPBTA_CNS_FIT_ONLY:-1}

# In CI we'll run an abbreviated version of the de novo signatures extraction
ABBREVIATED_MUTSIGS=${OPENPBTA_QUICK_MUTSIGS:-0}

# Split up the consensus MAF files by experimental strategy (writes to scratch)
Rscript --vanilla 02-split_experimental_strategy.R

if [ "$CNS_FIT_ONLY" == "0" ]; then

  # Run the mutational signatures analysis using existing signatures
  Rscript -e "rmarkdown::render('01-known_signatures.Rmd', clean = TRUE)"

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

fi

# Fit the 8 known CNS signatures
Rscript --vanilla 05-fit_cns_signatures.R  \
  --abbreviated $ABBREVIATED_MUTSIGS
  
# Process, visualize the 8 known CNS signature exposures
Rscript -e "rmarkdown::render('06-analyze_cns_fit.Rmd', clean = TRUE)"
