#!/bin/bash

# Josh Shapiro for CCDL 2019
#
# Runs 01-generate-independent-specimens.R with default settings.
# Takes one environment variable, `OPENPBTA_BASE_RELEASE`, if value is 1 then
# uses pbta-histologies-base.tsv, else runs all module with pbta-histologies.tsv (Default)

set -e
set -o pipefail

# If set to 1 (e.g., to generate files for a release), use the base histologies file
RUN_FOR_RELEASE=${OPENPBTA_BASE_RELEASE:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript -e "rmarkdown::render('00-repeated-samples.Rmd', params=list(base_run = ${RUN_FOR_RELEASE}), clean = TRUE)"

if [[ "$RUN_FOR_RELEASE" -eq "0" ]]
then
   HISTOLOGY_FILE="../../data/pbta-histologies.tsv" 
else 
   HISTOLOGY_FILE="../../data/pbta-histologies-base.tsv"  
fi

Rscript 01-generate-independent-specimens.R \
  -f $HISTOLOGY_FILE \
  -o results

# adding indepedent list of rnaseq samples
Rscript 02-generate-independent-rnaseq.R \
  --histology_file $HISTOLOGY_FILE \
  --output_directory results \
  --independent_dna_sample_df results/independent-specimens.wgswxs.primary-plus.tsv
