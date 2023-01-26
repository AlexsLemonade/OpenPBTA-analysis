#!/bin/bash

# Chante Bethell for CCDL 2019
#
# Run `01-filter-across-types.R`, `03-tumor-descriptor-and-assay-count.Rmd`, and '04-ns-for-manuscript.Rmd'
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

Rscript -e "rmarkdown::render('03-tumor-descriptor-and-assay-count.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('04-ns-for-manuscript.Rmd', clean = TRUE)"
