#!/bin/bash

# Pichai Raman - taken from --> (Chante Bethell for CCDL 2019)
#
# Run `01-run-ssgsea.R` and `02-analyze-plot-ssgsea-res.R` 
# sequentially. 

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript --vanilla 01-run-ssgsea.R
Rscript --vanilla 02-analyze-plot-ssgsea-res.R -a 1 -t 1 -p 1