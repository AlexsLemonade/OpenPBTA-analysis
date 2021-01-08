#!/bin/bash
# CCDL for ALSF 2021
# Joshua Shapiro
#
# Set so the script stops if there is an error
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# generate correlation matrix - rsem-tpm.polya
python3 01-correlation-matrix.py \
  ../../data/pbta-gene-expression-rsem-tpm.polya.rds \
  --clinical-path ../../data/pbta-histologies.tsv \
  --qc-manifest-path ../../data/pbta-mend-qc-manifest.tsv \
  --qc-results-path ../../data/pbta-mend-qc-results.tar.gz \
  --prefix rsem-tpm-polya- \
  --verbose

# generate correlation matrix - rsem-tpm.stranded
python3 01-correlation-matrix.py \
  ../../data/pbta-gene-expression-rsem-tpm.stranded.rds \
  --clinical-path ../../data/pbta-histologies.tsv \
  --qc-manifest-path ../../data/pbta-mend-qc-manifest.tsv \
  --qc-results-path ../../data/pbta-mend-qc-results.tar.gz \
  --prefix rsem-tpm-stranded- \
  --verbose

# generate thresholds and outliers - rsem-tpm.stranded
python3 02-thresholds-and-outliers.py \
  --prefix rsem-tpm-stranded- \
  --results results \
  --verbose
