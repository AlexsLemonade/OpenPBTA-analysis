####################################################
# Stephanie J. Spielman for ALSF CCDL 2020
#
# Run the GSEA pipeline, currently just `01-conduct-gsea-analysis.R`
# 
# Usage: bash run-gsea.sh
#
# Takes a single environmental variable `OPENPBTA_SMALLSET`. If FALSE or 0 (default), runs all samples. If TRUE or 1 runs a subset of samples (for testing in CI)
####################################################


set -e
set -o pipefail

SMALLSET=${OPENPBTA_SMALLSET:-1}


# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript --vanilla 01-conduct-gsea-analysis.R --smallset ${SMALLSET}
