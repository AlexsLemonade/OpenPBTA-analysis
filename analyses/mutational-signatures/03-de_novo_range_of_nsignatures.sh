#!/bin/bash

# JN Taroni for ALSF CCDL & SJ Spielman
# 2020

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# In CI we'll run an abbreviated version of the de novo signatures extraction
QUICK_MUTSIGS=${QUICK_MUTSIGS:-0}

# Which analysis are we running?
ANALYSIS=${ANALYSIS:-1}  # Expect 0 for GOF and 1 for Extraction


scratch_dir=../../scratch/mutational-signatures

# For initial GOF analysis to figure out how many k
gof_plot_dir=plots/denovo/gof
gof_result_dir=${scratch_dir}/gof

# For extraction once a limited range of k is assessed with GOF
extraction_plot_dir=plots/denovo/extraction
extraction_result_dir=${scratch_dir}/extraction


# Directories to hold all plots and results
mkdir -p $gof_plot_dir
mkdir -p $gof_result_dir
mkdir -p $extraction_plot_dir
mkdir -p $extraction_result_dir

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
  # GOF
  if [[ $ANALYSIS -eq "0" ]] 
  then
    FLOOR=2
    CEIL=8
    ITER=1000
    plot_dir=${gof_plot_dir}
    result_dir=${gof_result_dir}
  fi
  # Extraction
  if [[ $ANALYSIS -eq "1" ]] 
  then
    FLOOR=2
    CEIL=5
    ITER=3000
    plot_dir=${extraction_plot_dir}
    result_dir=${extraction_result_dir}
  fi  

  # Run sigfit with params
  for model in poisson multinomial; do
    for seed in {1..5}; do
       # De novo signatures extraction
       Rscript --vanilla \
         scripts/de_novo_signature_extraction.R \
         --maf_file ${maf_file} \
         --nsignatures_floor ${FLOOR} \
         --nsignatures_ceiling ${CEIL} \
         --num_iterations ${ITER} \
         --model ${model} \
         --seed ${seed} \
         --plot_output "${plot_dir}/seed_${seed}_model_${model}.png" \
         --output_file "${result_dir}/seed_${seed}_model_${model}.RDS"   
    done
  done
fi  
  
