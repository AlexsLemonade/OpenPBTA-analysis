#!/bin/bash

# Chante Bethell for CCDL 2020
#
# Run the HGG molecular subtyping pipeline.
# When OPENPBTA_SUBSET=1 (default), new HGG subset files will be generated.

set -e
set -o pipefail

# This option controls whether on not the step that generates the HGG only
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run the first script in this module that reclassifies high-grade gliomas
Rscript -e "rmarkdown::render('01-HGG-molecular-subtyping-defining-lesions.Rmd', clean = TRUE)"

# Run the second script in this module that subset files using the samples in the output
# file generated with `01-HGG-molecular-subtyping-defining-lesions.Rmd`.
if [ "$SUBSET" -gt "0" ]; then
  Rscript --vanilla 02-HGG-molecular-subtyping-subset-files.R
fi

# Run the copy number data cleaning notebook
Rscript -e "rmarkdown::render('03-HGG-molecular-subtyping-cnv.Rmd', clean = TRUE)"
