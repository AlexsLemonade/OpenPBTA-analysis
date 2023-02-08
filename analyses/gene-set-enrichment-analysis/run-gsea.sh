#########################################################################
# Stephanie J. Spielman for ALSF CCDL 2020
#
# Run the GSEA pipeline:
## 1. `01-conduct-gsea-analysis.R` to calculate scores
## 2. `02-model-gsea.Rmd` to model the scores
## 3. `03-assess-gsea-at-threshold.Rmd` to explore how results might be
#      affected by a tumor purity threshold on samples
#
# Usage: bash run-gsea.sh
#
# Takes two environment variables, `OPENPBTA_TESTING` and `OPENPBTA_BASE_SUBTYPING`.
#
# For `OPENPBTA_TESTING`:
#  If value is 1, the module is run with CI data
#  If value is 0, the module is run with the full OpenPBTA cohort (Default)
#
# For `OPENPBTA_BASE_SUBTYPING`:
#  If value is 1, only scripts required for subtyping are run
#  If value is 0, the full module is run (Default)
#########################################################################


set -e
set -o pipefail

IS_CI=${OPENPBTA_TESTING:-0}
RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

DATA_DIR="../../data"
RESULTS_DIR="results"

######## Calculate scores from polyA expression data ############
INPUT_FILE="${DATA_DIR}/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
OUTPUT_FILE="${RESULTS_DIR}/gsva_scores_polya.tsv"
Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE}

######## Calculate scores from stranded expression data ############
INPUT_FILE="${DATA_DIR}/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
OUTPUT_FILE="${RESULTS_DIR}/gsva_scores_stranded.tsv"
Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE}


# Only run the following if we are _not_ subtyping:
if [[ "$RUN_FOR_SUBTYPING" -lt "1" ]]; then
  ######## Model GSVA scores ############
  # Only run when pbta-histologies.tsv is generated which has harmonized_diagnosis
  Rscript -e "rmarkdown::render('02-model-gsea.Rmd', clean = TRUE, params=list(is_ci = ${IS_CI}))"

  ######## Calculate scores from stranded expression data for threshold-passing tumors only #######
  INPUT_FILE="${DATA_DIR}/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
  OUTPUT_FILE="${RESULTS_DIR}/gsva_scores_stranded_thresholded.tsv"
  Rscript --vanilla 01-conduct-gsea-analysis.R --input ${INPUT_FILE} --output ${OUTPUT_FILE} --apply_tumor_purity_threshold

  # Assess results generated for threshold-passing tumors
  Rscript -e "rmarkdown::render('03-assess-gsea-at-threshold.Rmd', clean = TRUE)"
fi










