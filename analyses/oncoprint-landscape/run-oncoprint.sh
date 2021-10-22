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

#### Files ---------------------------------------------------------------------

maf_consensus=../../data/pbta-snv-consensus-mutation.maf.tsv.gz
hotspots_maf=../../data/pbta-snv-scavenged-hotspots.maf.tsv.gz
fusion_file=../../data/pbta-fusion-putative-oncogenic.tsv
histologies_file=../../data/pbta-histologies.tsv
intermediate_directory=../../scratch/oncoprint_files
primary_filename="primary_only"
primaryplus_filename="primary-plus"
consensus_seg_autosomes_cnv_file=../../data/consensus_seg_annotated_cn_autosomes.tsv.gz
consensus_seg_cnv_xy_cnv_file=../../data/consensus_seg_annotated_cn_x_and_y.tsv.gz
oncoprint_data_directory=data

#### Prep genes of interest lists ----------------------------------------------

# Turn 1 manually curated CSV file that contains all histologies into
# individual TSV files
Rscript --vanilla 00-prepare-goi-lists.R

#### Map between DNA and RNA specimens -----------------------------------------

### Primary only samples mapping for oncoprint

Rscript --vanilla 01-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --hotspots_maf_file ${hotspots_maf} \
  --cnv_autosomes_file ${consensus_seg_autosomes_cnv_file} \
  --cnv_xy_file ${consensus_seg_cnv_xy_cnv_file} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primary_filename} \
  --independent_specimens ../../data/independent-specimens.wgs.primary.tsv

#### Primary plus samples mapping for oncoprint

Rscript --vanilla 01-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --hotspots_maf_file ${hotspots_maf} \
  --cnv_autosomes_file ${consensus_seg_autosomes_cnv_file} \
  --cnv_xy_file ${consensus_seg_cnv_xy_cnv_file} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primaryplus_filename} \
  --independent_specimens ../../data/independent-specimens.wgs.primary-plus.tsv

#### Oncoprints by broad histology ---------------------------------------------

# We'll use a declarative array to loop through pairs of broad histology labels
# and genes of interest files
histologies=(lgat embryonal hgat other)

declare -A labels=(
  [lgat]="Low-grade astrocytic tumor"
  [embryonal]="Embryonal tumor"
  [hgat]="Diffuse astrocytic and oligodendroglial tumor"
  [other]="Other CNS"
)

declare -A goi_files=(
  [lgat]="lgat_goi_list.tsv"
  [embryonal]="embryonal-tumor_goi_list.tsv"
  [hgat]="hgat_goi_list.tsv"
  [other]="other_goi_list.tsv"
)

# Will create two plots - primary only and "primary plus" samples
filenames=($primary_filename $primaryplus_filename)

# Print oncoprints by broad histology
for histology in "${histologies[@]}"; do
  for filename in "${filenames[@]}"; do

    # Print primary only oncoprints by broad histology
    Rscript --vanilla 02-plot-oncoprint.R \
      --maf_file "${intermediate_directory}/${filename}_maf.tsv" \
      --cnv_file "${intermediate_directory}/${filename}_cnv.tsv" \
      --fusion_file "${intermediate_directory}/${filename}_fusions.tsv" \
      --metadata_file "${histologies_file}" \
      --png_name "${filename}_${histology}_oncoprint.png" \
      --broad_histology "${labels[$histology]}" 

    # Genes of interest only version of oncoprint
    Rscript --vanilla 02-plot-oncoprint.R \
      --maf_file "${intermediate_directory}/${filename}_maf.tsv" \
      --cnv_file "${intermediate_directory}/${filename}_cnv.tsv" \
      --fusion_file "${intermediate_directory}/${filename}_fusions.tsv" \
      --metadata_file "${histologies_file}" \
      --goi_list "${oncoprint_data_directory}/${goi_files[$histology]}" \
      --top_n 20 \
      --png_name "${filename}_${histology}_goi_oncoprint.png" \
      --broad_histology "${labels[$histology]}" \
      --output_table "${filename}_${histology}_oncoprint_summary_n.tsv"

  done
done

Rscript --vanilla 03-oncoprint-n-count-table.R \
  --maf_file_po "${intermediate_directory}/primary_only_maf.tsv" \
  --cnv_file_po "${intermediate_directory}/primary_only_cnv.tsv" \
  --fusion_file_po "${intermediate_directory}/primary_only_fusions.tsv" \
  --maf_file_pp "${intermediate_directory}/primary-plus_maf.tsv" \
  --cnv_file_pp "${intermediate_directory}/primary-plus_cnv.tsv" \
  --fusion_file_pp "${intermediate_directory}/primary-plus_fusions.tsv" \
  --metadata_file "${histologies_file}" \
  --broad_histology_list "Low-grade astrocytic tumor,Embryonal tumor,Diffuse astrocytic and oligodendroglial tumor,Other CNS" \
  --short_name_list "lgat,embryonal,hgat,other" 
