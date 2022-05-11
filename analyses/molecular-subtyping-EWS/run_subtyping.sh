#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run notebook to subtype EWS per sample_id if  hallmark fusion in RNAseq samples 

Rscript -e "rmarkdown::render('01-run-subtyping-ewings.Rmd', clean = TRUE)"


