# C. Bethell and C. Savonen for CCDL 2019
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# Run original files - will not by default
RUN_ORIGINAL=${RUN_ORIGINAL:-0}

# Run testing files for circle CI - will not by default
IS_CI=${OPENPBTA_TESTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
results_dir=../../analyses/focal-cn-file-preparation/results
histologies_file=${data_dir}/histologies.tsv
gtf_file=${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz
goi_file=../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt
independent_specimens_file=${data_dir}/independent-specimens.wgswxs.primary.tsv

# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

# Run annotation step for consensus file
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${scratch_dir}/consensus_seg_with_status.tsv \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "consensus_seg_annotated_cn" \
--seg


# if we want to process the CNV data from the original callers
# (e.g., CNVkit, ControlFreeC)
if [ "$RUN_ORIGINAL" -gt "0" ]; then

# Prep the CNVkit data
Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Run annotation step for CNVkit
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${scratch_dir}/cnvkit_with_status.tsv \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "cnvkit_annotated_cn" \
--seg \
--runWXSonly

# Run annotation step for ControlFreeC
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${data_dir}/cnv-controlfreec.tsv.gz \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "controlfreec_annotated_cn" \
--controlfreec \
--runWXSonly

# filenameLead=("cnvkit_annotated_cn" "controlfreec_annotated_cn" "cnvkit_annotated_cn_wxs" "controlfreec_annotated_cn_wxs")
# chromosomeType=("autosomes" "x_and_y")
# for filename in ${filenameLead[@]}; do
#   for chromosome_type in ${chromosomesType[@]}; do
#       Rscript --vanilla rna-expression-validation.R \
#         --annotated_cnv_file results/${filename}_${chromosome_type}.tsv.gz \
#         --expression_file ${data_dir}/gene-expression-rsem-tpm-collapsed.rds \
#         --independent_specimens_file $independent_specimens_file \
#         --metadata $histologies_file \
#         --goi_list $goi_file \
#         --filename_lead ${filename}_${chromosome_type}
#   done
# done

fi