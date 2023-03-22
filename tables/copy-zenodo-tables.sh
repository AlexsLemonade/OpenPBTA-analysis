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

# Figure 3A and 3B: interaction plot

cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3a-data.csv" "${OUTDIR}/figure-3a-data.csv"
cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3b-data.csv" "${OUTDIR}/figure-3b-data.csv"

