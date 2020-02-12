#!/bin/bash

set -e
set -o pipefail

IS_CI=${OPENPBTA_TESTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript -e "rmarkdown::render('01-fusion-summary.Rmd', params=list(is_ci = ${IS_CI}), clean = TRUE)"
