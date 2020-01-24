#########################################################################
# Stephanie J. Spielman for ALSF CCDL 2020
#
# Run the GSEA pipeline:
## 1. `01-conduct-gsea-analysis.R` to calculate scores
## 2. `02-exploratory-gsea.Rmd` to explore the scores, lightly for now
# 
# Usage: bash run-gsea.sh 
#
# Takes one environment variable, `OPENPBTA_TESTING`, which is 1 for running
# samples in CI for testing, 0 for running the full dataset (Default)
#
#########################################################################


set -e
set -o pipefail


IS_CI=${OPENPBTA_TESTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


######## Calculate scores from polyA expression data ############
INPUT_FILE="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
OUTPUT_FILE="gsva_scores_polya.tsv"
Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE}


######## Calculate scores from stranded expression data ############
INPUT_FILE="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
OUTPUT_FILE="gsva_scores_stranded.tsv"
Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE}


######## Model GSVA scores ############
Rscript -e "rmarkdown::render('02-model-gsea.Rmd', clean = TRUE, params=list(is_ci = ${IS_CI}))" 
