#!/bin/bash

# Chante Bethell for CCDL 2019
#
# Run `01-ATRT-molecular-subtyping-data-prep.Rmd` and 
# `02-ATRT-molecular-subtyping-plotting.R` sequentially.

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript -e "rmarkdown::render('01-ATRT-molecular-subtyping-data-prep.Rmd', clean = TRUE)"
Rscript --vanilla 02-ATRT-molecular-subtyping-plotting.R
