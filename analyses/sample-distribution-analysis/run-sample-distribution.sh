#!/bin/bash

# Chante Bethell for CCDL 2019
#
# Run `01-filter-across-types.R` and `02-multilayer-plots.R`
# sequentially.

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript --vanilla 01-filter-across-types.R
Rscript --vanilla 02-multilayer-plots.R
