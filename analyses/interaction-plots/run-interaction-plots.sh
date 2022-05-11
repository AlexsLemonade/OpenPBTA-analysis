#########################################################################
# Module author: JA Shapiro for ALSF CCDL 2019-2020
# Script author: Jo Lynne Rokita for D3b 2022
#
# Run the interaction plots pipeline:
## 1. `01-create-interaction-plots` to create interaction plots
## 2. `02-result-ns-for-manuscript.Rmd` to create sample Ns used for the manuscript
#
# Usage: bash run-interaction-plots.sh
#########################################################################

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run create interaction plots
bash 01-create-interaction-plots.sh

# summarize output from both classifiers and expected classification
Rscript -e "rmarkdown::render('02-result-ns-for-manuscript.Rmd', clean = TRUE)"