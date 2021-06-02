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

### Collapse rsem-expected_counts using 
### OpenPBTA-analysis/analyses/collapse-rnaseq module

# create results directory if it doesn't already exist
mkdir -p results
mkdir -p plots

echo 'Run TMM normalized differential gene expression and stably expressed gene analysis...'
Rscript --vanilla '01-protocol-dge-seg.R' -n 'tmm'

echo 'Run UQ-pgQ2 normalized differential gene expression and stably expressed gene analysis...'
Rscript --vanilla '01-protocol-dge-seg.R' -n 'uqpgq2'

echo 'Run DESeq2 standard differential gene expression analysis...'
Rscript --vanilla '02-deseq2-protocol-dge.R'
