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
denovo_plot_dir=plots/denovo/gof
denovo_results_dir=${scratch_dir}/gof

# Directories to hold the goodness-of-fit plots and results
mkdir -p $denovo_plot_dir
mkdir -p $denovo_results_dir

# The MAF file we'll use is going to WGS samples only
maf_file=${scratch_dir}/pbta-snv-consensus-wgs.tsv.gz

# If abbreviated mutsigs is "true", we'll extract 3 signatures using
# an unacceptably low number of iterations to do any analysis with
if [ "$QUICK_MUTSIGS" -gt "0" ]
then
  # De novo signatures extraction
  Rscript --vanilla \
    scripts/de_novo_signature_extraction.R \
    --maf_file ${maf_file} \
    --nsignatures_floor 3 \
    --nsignatures_ceiling 4 \
    --num_iterations 10 \
    --seed 42 
else
  for model in multinomial poisson; do
    for seed in {1..5}; do

       # De novo signatures extraction
       Rscript --vanilla \
         scripts/de_novo_signature_extraction.R \
         --maf_file ${maf_file} \
         --nsignatures_floor 2 \
         --nsignatures_ceiling 8 \
         --num_iterations 1000 \
         --model ${model} \
         --seed ${seed} \
         --plot_output "${denovo_plot_dir}/gof_seed_${seed}_model_${model}.pdf" \
         --output_file "${denovo_results_dir}/gof_seed_${seed}_model_${model}.RDS"
      

    done
  done

fi
