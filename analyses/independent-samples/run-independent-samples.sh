#!/bin/bash

# Josh Shapiro for CCDL 2019
#
# Runs 01-generate-independent-specimens.R with default settings.
# Takes one environment variable, `BASE_SUBTYPING`, if value is 1 then
# uses pbta-histologies-base.tsv for subtyping if value is 0 runs all modules with pbta-histologies.tsv(Default)

set -e
set -o pipefail

RUN_FOR_SUBTYPING=${BASE_SUBTYPING:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript -e "rmarkdown::render('00-repeated-samples.Rmd',params=list(base_run = ${RUN_FOR_SUBTYPING}), clean = TRUE)"

if [[ RUN_FOR_SUBTYPING == "0" ]]
then
Rscript 01-generate-independent-specimens.R \
  -f ../../data/pbta-histologies.tsv \
  -o results

# adding indepedent list of rnaseq samples
Rscript 02-generate-independent-rnaseq.R \
  --histology_file ../../data/pbta-histologies.tsv \
  --output_directory results \
  --independent_dna_sample_df ../../data/independent-specimens.wgswxs.primary-plus.tsv  

else 
Rscript 01-generate-independent-specimens.R \
  -f ../../data/pbta-histologies-base.tsv \
  -o results

# adding indepedent list of rnaseq samples
Rscript 02-generate-independent-rnaseq.R \
  --histology_file ../../data/pbta-histologies-base.tsv \
  --output_directory results \
  --independent_dna_sample_df ../../data/independent-specimens.wgswxs.primary-plus.tsv
fi
