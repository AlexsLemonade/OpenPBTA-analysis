#!/bin/bash

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir=../../data

# Run the mutational signatures analysis using existing signatures
Rscript -e "rmarkdown::render('01-mutational_signatures.Rmd', clean = TRUE)"

# TODO: this probably needs specific logic re: CI -- can drop number of iterations?
# Run the de novo mutational signatures step for a range of nsignatures values
Rscript --vanilla \
  scripts/de_novo_signature_fitting.R \
  --maf_file ${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz \
  --nsignatures_floor 5 \
  --nsignatures_ceiling 15 \
  --num_iterations 1000 \
  --seed 42 \
  --output_file results/sigfit_signatures.RDS
