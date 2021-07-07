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

Rscript -e "rmarkdown::render('00-repeated-samples.Rmd',params=list(base_run = ${RUN_FOR_SUBTYPING}), clean = TRUE)"

if [[ RUN_FOR_SUBTYPING -eq "0" ]]
then
   HISTOLOGY_FILE="../../data/histologies.tsv" 
else 
   HISTOLOGY_FILE="../../data/histologies-base.tsv"  
fi

independent_DNA_each="results/independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv"
independent_DNA_all="results/independent-specimens.wgswxspanel.primary-plus.allcohorts.tsv"

Rscript 01-generate-independent-specimens.R \
  -f $HISTOLOGY_FILE \
  -o results

# adding indepedent list of rnaseq samples
Rscript 02-generate-independent-rnaseq.R \
  --histology_file $HISTOLOGY_FILE \
  --output_directory results \
  --independent_dna_sample_df_each $independent_DNA_each \
  --independent_dna_sample_df_all $independent_DNA_all
