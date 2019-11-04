#!/bin/bash

# JA Shapiro for CCDL 2019
#
# Runs scripts/01-process_mutations.R with some default settings.

set -e
set -o pipefail

base_dir=analyses/interaction-plots
script_dir=${base_dir}/scripts
results_dir=${base_dir}/results
plot_dir=${base_dir}/plots


ind_samples=data/independent-specimens.wgs.primary-plus.tsv

# make output directories if they don't exist
mkdir -p $results_dir
mkdir -p $plot_dir

# run scripts
# using lancet data for now
lancet_cooccur=${results_dir}/lancet_top50.tsv
lancet_plot=${plot_dir}/lancet_top50.png

Rscript ${script_dir}/01-process_mutations.R \
  --maf data/pbta-snv-lancet.vep.maf.gz \
  --metadata data/pbta-histologies.tsv \
  --specimen_list ${ind_samples} \
  --vaf 0.2 \
  --out ${lancet_cooccur}

Rscript ${script_dir}/02-plot_interactions.R \
  --infile ${lancet_cooccur} \
  --outfile ${lancet_plot}