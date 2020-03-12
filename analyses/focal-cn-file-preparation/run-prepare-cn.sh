# C. Bethell and C. Savonen for CCDL 2019
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# Run original files - will not by default
RUN_ORIGINAL=${RUN_ORIGINAL:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
results_dir=../../analyses/focal-cn-file-preparation/results
histologies_file=${data_dir}/pbta-histologies.tsv
gtf_file=${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz
goi_file=../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt
independent_specimens_file=${data_dir}/independent-specimens.wgswxs.primary.tsv
ucsc_bed_file=${results_dir}/ucsc_cytoband.bed
consensus_bed_file=${scratch_dir}/consensus_seg_with_status.bed
loss_intersect_with_cytoband_file=${scratch_dir}/intersect_with_cytoband_losses.bed
gain_intersect_with_cytoband_file=${scratch_dir}/intersect_with_cytoband_gains.bed
callable_intersect_with_cytoband_file=${scratch_dir}/intersect_with_cytoband_callable.bed

# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

# Download and save UCSC cytoband file as bed file
wget -O ${scratch_dir}/ucsc_cytoband.bed http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
    
# Use bedtools intersect to find the intersection between the UCSC file with
# cytoband data and the `scratch/consensus_with_status.tsv` file prepared in
# `02-add-ploidy-consensus.Rmd`

bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b ${scratch_dir}/consensus_seg_with_status_losses.bed \
    -sorted \
    > $loss_intersect_with_cytoband_file

bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b ${scratch_dir}/consensus_seg_with_status_gains.bed \
    -sorted \
    > $gain_intersect_with_cytoband_file

bedtools coverage \
    -a ${scratch_dir}/ucsc_cytoband.bed \
    -b $consensus_bed_file \
    -sorted \
    > $callable_intersect_with_cytoband_file

# # Run annotation step for consensus file
# Rscript --vanilla 03-prepare-cn-file.R \
#   --cnv_file ${scratch_dir}/consensus_seg_with_status.tsv \
#   --gtf_file $gtf_file \
#   --metadata $histologies_file \
#   --filename_lead "consensus_seg_annotated_cn" \
#   --seg
# 
# libraryStrategies=("polya" "stranded")
# chromosomesType=("autosomes" "x_and_y")
# for strategy in ${libraryStrategies[@]}; do
# 
#   for chromosome_type in ${chromosomesType[@]}; do
# 
#     Rscript --vanilla rna-expression-validation.R \
#       --annotated_cnv_file results/consensus_seg_annotated_cn_${chromosome_type}.tsv.gz \
#       --expression_file ${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.${strategy}.rds \
#       --independent_specimens_file $independent_specimens_file \
#       --metadata $histologies_file \
#       --goi_list $goi_file \
#       --filename_lead "consensus_seg_annotated_cn"_${chromosome_type}_${strategy}
#   done
# done
# 
# # if we want to process the CNV data from the original callers
# # (e.g., CNVkit, ControlFreeC)
# if [ "$RUN_ORIGINAL" -gt "0" ]; then
# 
#   # Prep the CNVkit data
#   Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"
# 
#   # Run annotation step for CNVkit
#   Rscript --vanilla 03-prepare-cn-file.R \
#     --cnv_file ${scratch_dir}/cnvkit_with_status.tsv \
#     --gtf_file $gtf_file \
#     --metadata $histologies_file \
#     --filename_lead "cnvkit_annotated_cn" \
#     --seg
# 
#   # Run annotation step for ControlFreeC
#   Rscript --vanilla 03-prepare-cn-file.R \
#     --cnv_file ${data_dir}/pbta-cnv-controlfreec.tsv.gz \
#     --gtf_file $gtf_file \
#     --metadata $histologies_file \
#     --filename_lead "controlfreec_annotated_cn" \
#     --controlfreec
# 
#   filenameLead=("cnvkit_annotated_cn" "controlfreec_annotated_cn")
#   for filename in ${filenameLead[@]}; do
#     for strategy in ${libraryStrategies[@]}; do
#       for chromosome_type in ${chromosomesType[@]}; do
#         Rscript --vanilla rna-expression-validation.R \
#           --annotated_cnv_file results/${filename}_${chromosome_type}.tsv.gz \
#           --expression_file ${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.${strategy}.rds \
#           --independent_specimens_file $independent_specimens_file \
#           --metadata $histologies_file \
#           --goi_list $goi_file \
#           --filename_lead ${filename}_${chromosome_type}_${strategy}
#       done
#     done
#   done
# 
# fi
