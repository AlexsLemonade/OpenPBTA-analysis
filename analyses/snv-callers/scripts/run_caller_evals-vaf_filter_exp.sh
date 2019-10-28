#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose:Run TMB calculations at various VAF cutoffs for vaf_cutoff_experiment.Rmd

cd ../.. 

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# VAF cutoffs to use
vaf_cutoff=("0.10" "0.20" "0.30" "0.40")

# Same files as originally run
datasets=("strelka2" "mutect2" "lancet" "vardict")

# WGS files for each caller
wgs_files=("WGS.hg38.strelka2.unpadded.bed" "WGS.hg38.mutect2.unpadded.bed" "WGS.hg38.lancet.300bp_padded.bed" "WGS.hg38.vardict.100bp_padded.bed")

for ((i=0;i<${#vaf_cutoff[@]};i++)); 
  do
  for ((j=0;j<${#datasets[@]};j++)); 
    do
    echo "Processing dataset: ${datasets[$j]}"
    Rscript analyses/snv-callers/scripts/01-calculate_vaf_tmb.R \
      --label ${datasets[$j]} \
      --output analyses/snv-callers/results/vaf_filter/cutoff_${vaf_cutoff[$i]} \
      --file_format rds \
      --maf data/pbta-snv-${datasets[$j]}.vep.maf.gz \
      --metadata data/pbta-histologies.tsv \
      --bed_wgs data/${wgs_files[$j]} \
      --bed_wxs data/WXS.hg38.100bp_padded.bed \
      --vaf_filter ${vaf_cutoff[$i]} \
      --overwrite \
      --no_region
    done
  done
