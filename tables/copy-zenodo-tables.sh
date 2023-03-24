#!/bin/bash
#
# S. Spielman for CCDL, 2023
# This script copies figure data CSV files that were exported
#  within individual modules into the `zenodo-upload/` directory here.
#  Context: https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1692


# This script must be called from `run-tables.sh`, or it must be run from this directory.

# enviroment settings
set -eo pipefail

# Define input, output directories
ANALYSIS_DIR="../analyses"

OUTDIR="zenodo-upload"
mkdir -p "${OUTDIR}" # make directory in case

######################### Copy CSV files #########################

# Figure 3A and 3B: interaction and cooccurrence plots
cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3a-data.csv" "${OUTDIR}/figure-3a-data.csv"
cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3b-data.csv" "${OUTDIR}/figure-3b-data.csv"

# Figure 3C: Chromothripsis scatterplot
cp "${ANALYSIS_DIR}/chromothripsis/results/figure-3c-data.csv" "${OUTDIR}/figure-3c-data.csv"

# Figure S3A: Violin plot of tumor purity across cancer groups
cp "${ANALYSIS_DIR}/tumor-purity-exploration/results/figure-S3a-data.csv" "${OUTDIR}/figure-S3a-data.csv"


# Figure S7C: Sequencing center UMAP
cp "${ANALYSIS_DIR}/transcriptomic-dimension-reduction/results/figure-S7c-data.csv" "${OUTDIR}/figure-S7c-data.csv"

# Figures S7D, S7E, S7F: TP53 tumor purity thresholded re-analyses for TP53 results
cp "${ANALYSIS_DIR}/tp53_nf1_score/results/tumor-purity-threshold/figure-S7d-data.csv" "${OUTDIR}/figure-S7d-data.csv"
cp "${ANALYSIS_DIR}/tp53_nf1_score/results/tumor-purity-threshold/figure-S7e-data.csv" "${OUTDIR}/figure-S7e-data.csv"
cp "${ANALYSIS_DIR}/tp53_nf1_score/results/tumor-purity-threshold/figure-S7f-data.csv" "${OUTDIR}/figure-S7f-data.csv"

# Figure S7H: Tumor purity thresholded UMAP
cp "${ANALYSIS_DIR}/transcriptomic-dimension-reduction/results/figure-S7h-data.csv" "${OUTDIR}/figure-S7h-data.csv"

# Figure S7I: Tumor purity thresholded quanTIseq results
cp "${ANALYSIS_DIR}/immune-deconv/results/figure-S7i-data.csv" "${OUTDIR}/figure-S7i-data.csv"
