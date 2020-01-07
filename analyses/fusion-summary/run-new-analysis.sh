#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

data_dir="../../data"
scratch_dir="../..scratch"

demographic_file="${data_dir}/pbta-histologies.tsv"
fusions_file="${data_dir}/pbta-fusion-putative-oncogenic.tsv"

Rscript --vanilla ${script_directory}/01-fusion-summary.R \
  --demographic_file ${demographic_file} \
  --fusions_file ${fusions_file} \
  --output_dir ${script_directory}/results 
