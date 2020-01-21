#!/bin/bash
#
# Jaclyn Taroni for ALSF CCDL 2020
#
# This shell script runs the analysis module for molecularly subtyping embryonal
# tumors.

set -e
set -o pipefail

# This option controls whether on not the step that generates the subset
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run the first script in this module that identifies non-ATRT and non-MB
# embryonal tumors and those tumors with TTYH1 fusions for the purposes of
# subsetting files downstream
Rscript -e "rmarkdown::render('01-samples-to-subset.Rmd', clean = TRUE)"

# Run the second script in this module that subset files using the samples in
# the output file generated with `01-samples-to-subset.Rmd`.
if [ "$SUBSET" -gt "0" ]; then
  Rscript --vanilla 02-generate-subset-files.R
fi
