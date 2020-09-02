#!/bin/bash
# C. Savonen
# CCDL for ALSF 2020

# Purpose: Calculate TMB without taking the nonsynonymous filter for TCGA data
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# The sqlite database made from the callers will be called:
pbta_dbfile=scratch/snv_db.sqlite
tcga_dbfile=scratch/tcga_snv_db.sqlite

# BED and GTF file paths
cds_file=scratch/gencode.v27.primary_assembly.annotation.bed

# Calculate consensus TMB without nonsynfilter

# For PBTA
Rscript ../scripts/03-calculate_tmb.R \
  --db_file $pbta_dbfile \
  --output analyses/snv-callers/results/no_filter \
  --metadata data/pbta-histologies.tsv \
  --coding_regions $cds_file \
  --overwrite

# For TCGA
Rscript ../scripts/03-calculate_tmb.R \
  --db_file $tcga_dbfile \
  --output analyses/snv-callers/results/no_filter \
  --metadata data/pbta-tcga-manifest.tsv \
  --coding_regions $cds_file \
  --overwrite \
  --tcga

# Run the notebook that makes plots
Rscript -e "rmarkdown::render('explore_nonsynfilter.Rmd', clean = TRUE)"
