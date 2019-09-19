# Bethell and Taroni for CCDL 2019
# Generates lists of scatter plots for dimension reduction techniques.
# 
# Usage: bash 03-get-dimension-reduction-plot-lists.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

INPUT="results"
OUTPUT="plots"

#### Broad histology plots -----------------------------------------------------

declare -a arr=("rsem_all" "rsem_stranded" "rsem_polyA" "kallisto_all" "kallisto_stranded" "kallisto_polyA")

for filename_lead in "${arr[@]}"
do
  Rscript --vanilla scripts/get-plot-list.R \
    --input_directory ${INPUT} \
    --filename_lead ${filename_lead} \
    --output_directory ${OUTPUT} \
    --color_variable broad_histology
done

#### Library strategy plots ----------------------------------------------------

Rscript --vanilla scripts/get-plot-list.R \
  --input_directory ${INPUT} \
  --filename_lead rsem_all \
  --output_directory ${OUTPUT} \
  --color_variable RNA_library

Rscript --vanilla scripts/get-plot-list.R \
  --input_directory ${INPUT} \
  --filename_lead kallisto_all \
  --output_directory ${OUTPUT} \
  --color_variable RNA_library
