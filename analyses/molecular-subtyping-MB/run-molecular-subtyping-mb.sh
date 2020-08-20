#!/bin/bash
# Module author: Komal S. Rathi
# 2020

# This script runs the steps for molecular subtyping of Medulloblastoma samples

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# filter to MB samples and/or batch correct
Rscript --vanilla 00-filter-and-batch-correction.R \
--batch_col RNA_library \
--output_prefix medulloblastoma-exprs \
--output_dir input

# classify MB subtypes
Rscript --vanilla 01-classify-mb.R \
--corrected_mat input/medulloblastoma-exprs-batch-corrected.rds \
--uncorrected_mat input/medulloblastoma-exprs.rds \
--output_prefix mb-classified

# summarize expected and observed classification
Rscript -e "rmarkdown::render('02-compare-classes.Rmd', clean = TRUE)"
