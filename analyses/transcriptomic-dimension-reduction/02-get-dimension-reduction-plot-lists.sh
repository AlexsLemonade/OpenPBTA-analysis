# Bethell and Taroni for CCDL 2019
# Generates lists of scatter plots for dimension reduction techniques.
#
# Usage: bash 02-get-dimension-reduction-plot-lists.sh

COLORVAR=${COLOR:-broad_histology}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

INPUT="results"
OUTPUT="plots/plot_data"

#### Broad histology plots -----------------------------------------------------

declare -a arr=("rsem_stranded_none" "rsem_stranded_log" "rsem_polyA_none" "rsem_polyA_log" "kallisto_stranded_none" "kallisto_stranded_log" "kallisto_polyA_none" "kallisto_polyA_log")

for filename_lead in "${arr[@]}"
do
  Rscript --vanilla scripts/get-plot-list.R \
    --input_directory ${INPUT} \
    --filename_lead ${filename_lead} \
    --output_directory ${OUTPUT} \
    --color_variable ${COLORVAR}
done
