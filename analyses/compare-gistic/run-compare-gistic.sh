#!/bin/bash

# Chante Bethell for CCDL 2020
#
# Run the GISTIC comparison R notebooks in this module sequentially.

set -e
set -o pipefail

IS_CI=${OPENPBTA_TESTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript -e "rmarkdown::render('01-GISTIC-cohort-vs-histology-comparison.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('02-GISTIC-tidy-data-prep.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('03-GISTIC-gene-level-tally.Rmd', clean = TRUE, params = list(is_ci = ${IS_CI}))"
