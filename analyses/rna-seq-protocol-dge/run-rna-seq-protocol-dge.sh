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

# generate collapsed matrices for poly-A and stranded datasets
libraryStrategies=("polya" "stranded")
for strategy in ${libraryStrategies[@]}; do

  Rscript --vanilla 01-summarize_matrices.R \
    -i ../../data/pbta-gene-counts-rsem-expected_count.${strategy}.rds \
    -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
    -m results/pbta-gene-counts-rsem-expected_count-collapsed.${strategy}.rds \
    -t results/pbta-gene-counts-rsem-expected_count-collapsed_table.${strategy}.rds

done
