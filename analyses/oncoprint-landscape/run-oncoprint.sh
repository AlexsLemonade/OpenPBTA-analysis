# Chante Bethell for CCDL 2019
# Run 01-plot-oncoprint.R
#
# Usage: bash run-oncoprint.sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# For the genes lists
# https://stackoverflow.com/questions/1527049/how-can-i-join-elements-of-an-array-in-bash
function join_by { local IFS="$1"; shift; echo "$*"; }

#### Files

maf_consensus=../../data/pbta-snv-consensus-mutation.maf.tsv.gz
fusion_file=../../data/pbta-fusion-putative-oncogenic.tsv
histologies_file=../../data/pbta-histologies.tsv
intermediate_directory=../../scratch/oncoprint_files
primary_filename="primary_only"
primaryplus_filename="primary-plus"
focal_directory=../focal-cn-file-preparation/results
focal_cnv_file=${focal_directory}/consensus_seg_most_focal_cn_status.tsv.gz

# each element of the array is a file that contains genes of interest
genes_list=("../interaction-plots/results/gene_disease_top50.tsv" \
            "../focal-cn-file-preparation/results/consensus_seg_focal_cn_recurrent_genes.tsv")
# join into a string, where file paths are separated by commas
genes_list=$(join_by , "${genes_list[@]}")

### Primary only samples mapping for oncoprint

Rscript --vanilla 00-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --cnv_file ${focal_cnv_file} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primary_filename} \
  --independent_specimens ../../data/independent-specimens.wgs.primary.tsv

#### Primary plus samples mapping for oncoprint

Rscript --vanilla 00-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --cnv_file ${focal_cnv_file} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primaryplus_filename} \
  --independent_specimens ../../data/independent-specimens.wgs.primary-plus.tsv

# Print oncoprints by broad histology
for histology in "Low-grade astrocytic tumor" \
"Embryonal tumor" \
"Diffuse astrocytic and oligodendroglial tumor" \
"Ependymal tumor" \
"Other CNS"
do
  
  # Print primary only oncoprints by broad histology
  Rscript --vanilla 01-plot-oncoprint.R \
    --maf_file ${intermediate_directory}/${primary_filename}_maf.tsv \
    --cnv_file ${intermediate_directory}/${primary_filename}_cnv.tsv \
    --fusion_file ${intermediate_directory}/${primary_filename}_fusions.tsv \
    --metadata_file ${histologies_file} \
    --png_name ${primary_filename}_"$histology"_oncoprint.png \
    --broad_histology "$histology"
      
  # Genes of interest only version of oncoprint
  Rscript --vanilla 01-plot-oncoprint.R \
    --maf_file ${intermediate_directory}/${primary_filename}_maf.tsv \
    --cnv_file ${intermediate_directory}/${primary_filename}_cnv.tsv \
    --fusion_file ${intermediate_directory}/${primary_filename}_fusions.tsv \
    --metadata_file ${histologies_file} \
    --goi_list ${genes_list} \
    --png_name ${primary_filename}_"${histology}"_goi_oncoprint.png \
    --broad_histology "${histology}"


  # Print primary plus oncoprints by broad histology
  Rscript --vanilla 01-plot-oncoprint.R \
    --maf_file ${intermediate_directory}/${primaryplus_filename}_maf.tsv \
    --cnv_file ${intermediate_directory}/${primaryplus_filename}_cnv.tsv \
    --fusion_file ${intermediate_directory}/${primaryplus_filename}_fusions.tsv \
    --metadata_file ${histologies_file} \
    --png_name ${primaryplus_filename}_"$histology"_oncoprint.png \
    --broad_histology "$histology"

  # Genes of interest only version of oncoprint
  Rscript --vanilla 01-plot-oncoprint.R \
    --maf_file ${intermediate_directory}/${primaryplus_filename}_maf.tsv \
    --cnv_file ${intermediate_directory}/${primaryplus_filename}_cnv.tsv \
    --fusion_file ${intermediate_directory}/${primaryplus_filename}_fusions.tsv \
    --metadata_file ${histologies_file} \
    --goi_list ${genes_list} \
    --png_name ${primaryplus_filename}_"$histology"_goi_oncoprint.png \
    --broad_histology "$histology"
    
done
