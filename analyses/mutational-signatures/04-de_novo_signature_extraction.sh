#!/bin/bash

# JN Taroni for ALSF CCDL & SJ Spielman
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# In CI we'll run an abbreviated version of the de novo signatures extraction
QUICK_MUTSIGS=${QUICK_MUTSIGS:-0}

scratch_dir=../../scratch/mutational-signatures
denovo_results_dir=${scratch_dir}/extraction

# Directories to hold the results
mkdir -p $denovo_results_dir

# The MAF file we'll use is going to WGS samples only
maf_file=${scratch_dir}/pbta-snv-consensus-wgs.tsv.gz

# If abbreviated mutsigs is "true", we'll use an 
# unacceptably low number of iterations to do any analysis
if [ "$QUICK_MUTSIGS" -gt "0" ]
then
    Rscript --vanilla \
      scripts/de_novo_signature_extraction.R \
      --maf_file ${maf_file} \
      --nsignatures_floor 3 \
      --nsignatures_ceiling 3 \
      --num_iterations 10 \
      --seed 42 
else 
    # De novo signatures extraction with \k in 3:5 for each of two models
    for model in multinomial poisson; do
        for K in {3..5}; do
            Rscript --vanilla \
              scripts/de_novo_signature_extraction.R \
              --maf_file ${maf_file} \
              --nsignatures_floor ${K} \
              --nsignatures_ceiling ${K} \
              --num_iterations 3000 \
              --model ${model} \
              --seed 42 \
              --output_file "${denovo_results_dir}/denovo_k_${K}_model_${model}.RDS"
    done
  done
fi
