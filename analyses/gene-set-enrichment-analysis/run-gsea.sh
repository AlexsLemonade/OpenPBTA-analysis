#########################################################################
# Stephanie J. Spielman for ALSF CCDL 2020
#
# Run the GSEA pipeline, currently just `01-conduct-gsea-analysis.R`
# 
# Usage: bash run-gsea.sh
#
#########################################################################


set -e
set -o pipefail


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
