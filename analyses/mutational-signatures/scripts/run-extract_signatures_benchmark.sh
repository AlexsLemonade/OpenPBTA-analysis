#!/bin/bash

# Stephanie Spielman, 2021 
# This script does not need to be run in CI. It will take a long time (few days) and requires about 20 GB RAM.

set -e
set -o pipefail

# Set the working directory to the directory of this file
WDIR="$(dirname "${BASH_SOURCE[0]}")"
cd $WDIR

# Create output path for goodness-of-fit plots, if not created
PLOT_DIR=../plots/denovo/gof
mkdir -p ${PLOT_DIR}

# Grab WGS maf file from scratch (produced by 02-split_experimental_strategy.R)
MAF_FILE=${WDIR}/../../../scratch/mutational-signatures/pbta-snv-consensus-wgs.tsv.gz

# Run benchmarking across 5 seeds and 2 models
for model in multinomial poisson; do
  for seed in {1..5}; do
    Rscript --vanilla \
      extract_signatures_benchmark.R \
      --maf_file ${MAF_FILE} \
      --model ${model} \
      --seed ${seed} \
      --plot_dir ${PLOT_DIR} 
  done
done
