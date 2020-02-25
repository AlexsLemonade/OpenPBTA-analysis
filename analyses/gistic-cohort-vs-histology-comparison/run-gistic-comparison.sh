#!/bin/bash

# Chante Bethell for CCDL 2020
#
# Run `gistic-cohort-vs-histology-comparison.Rmd` and
# `gistic-gene-level-call-comparison.Rmd` sequentially.

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript -e "rmarkdown::render('gistic-cohort-vs-histology-comparison.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('gistic-gene-level-call-comparison.Rmd', clean = TRUE)"