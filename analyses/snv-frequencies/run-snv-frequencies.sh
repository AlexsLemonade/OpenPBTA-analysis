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

Rscript --vanilla '01-snv-frequencies.R'

# --no-name option stops the filename and timestamp from being stored in the
# output file. So rerun will have the same file.
gzip --no-name results/snv-consensus-annotated-mut-freq.json
