#!/bin/bash
# PediatricOpenTargets 2021
# Yuanchao Zhang
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# copied from https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/scripts/run_in_ci.sh
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# create results directory if it doesn't already exist
mkdir -p results
mkdir -p plots

echo 'Prepare data...'
Rscript --vanilla '00-prepare-data.R'

echo 'Run RUVSeq DESeq2 differential gene expression analysis on RNA-seq libraries with matched sample IDs...'
Rscript --vanilla '01-protocol-ruvseq.R' -d 'match'

echo 'Run RUVSeq DESeq2 differential gene expression analysis on DIPG RNA-seq libraries without matching sample IDs...'
Rscript --vanilla '01-protocol-ruvseq.R' -d 'dipg'

echo 'Run RUVSeq DESeq2 differential gene expression analysis on NBL RNA-seq libraries...'
Rscript --vanilla '01-protocol-ruvseq.R' -d 'nbl'
