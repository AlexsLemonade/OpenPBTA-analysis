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

# If results directory exists, remove the existing results.
# If existing .gz files exist, gzip command asks for confirmation.
# Design-wise, it is not surprising to remove previously generated results
# before a new run.
if [[ -d results ]]; then
    rm -r results
fi

mkdir -p results

Rscript --vanilla '01-tpm-summary-stats.R'

# The --no-name option stops the filename and timestamp from being stored in the
# output file. So rerun will have the same file.
gzip --no-name results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv
gzip --no-name results/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv
gzip --no-name results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.json
gzip --no-name results/long_n_tpm_mean_sd_quantile_group_wise_zscore.json
