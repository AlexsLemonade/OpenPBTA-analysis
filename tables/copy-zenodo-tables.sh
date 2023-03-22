#!/bin/bash
#
# S. Spielman for CCDL, 2023
# This script copies figure data CSV files that were exported
#  within individual modules into the `zenodo-upload/` directory here.
#  Context: https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1692

# enviroment settings
set -eo pipefail

# Find directory of this script
WORKDIR="$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
BASEDIR="$(dirname "$WORKDIR")"

# Define input, output directories
ANALYSIS_DIR="${BASEDIR}/analyses"

OUTDIR="${BASEDIR}/tables/zenodo-upload"
mkdir -p "${OUTDIR}" # make directory in case



######################### Copy CSV files #########################

# Figure 3A and 3B: interaction plot

cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3a-data.csv" "${OUTDIR}/figure-3a-data.csv"
cp "${ANALYSIS_DIR}/interaction-plots/results/figure-3b-data.csv" "${OUTDIR}/figure-3b-data.csv"


