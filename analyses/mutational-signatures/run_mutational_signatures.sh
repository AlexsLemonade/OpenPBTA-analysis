#!/bin/bash

# J. Taroni for ALSF CCDL 
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# In CI we'll run an abbreviated version of the de novo signatures extraction
ABBREVIATED_MUTSIGS=${OPENPBTA_QUICK_MUTSIGS:-0}

data_dir=../../data

# Run the mutational signatures analysis using existing signatures
Rscript -e "rmarkdown::render('01-mutational_signatures.Rmd', clean = TRUE)"

# if abbreviated mutsigs is "true", we'll run only one number of signatures and
# for an unacceptably low number of iterations to do any analysis with
if [ "$ABBREVIATED_MUTSIGS" -gt "0" ]
then
  FLOOR=10
  CEILING=10
  NUM_ITER=10
else
  FLOOR=5
  CEILING=15
  NUM_ITER=1000
fi

# De novo signatures extraction
Rscript --vanilla \
  scripts/de_novo_signature_fitting.R \
  --maf_file ${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz \
  --nsignatures_floor $FLOOR \
  --nsignatures_ceiling $CEILING \
  --num_iterations $NUM_ITER \
  --seed 42 \
  --output_file results/denovo_sigfit_signatures.RDS

# Notebook for picking k and deriving k = 10 files - will not get run in CI
Rscript -e "rmarkdown::render('02-explore_de_novo.Rmd', clean = TRUE, params = list(is_ci = ${ABBREVIATED_MUTSIGS}))"

# Notebook for matching de novo signatures to reference signatures
Rscript -e "rmarkdown::render('03-match_de_novo.Rmd', clean = TRUE)"

# Notebook for de novo signatures per sample and plotting
Rscript -e "rmarkdown::render('04-de_novo_per_sample.Rmd', clean = TRUE)"
