#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: 1) Run an intial evaluation of each variant caller's MAF file
#          2) Compare the variant callers to each other. (Coming soon)

# Change directory
cd kitematic

# The files named in these arrays will be ran in the analysis. 
declare -a datasets=("strelka2" "mutect2" "lancet" "vardict")
declare -a wgs_files=("WGS.hg38.strelka2.unpadded.bed" "WGS.hg38.mutect2.unpadded.bed" "WGS.hg38.lancet.unpadded.bed" "WGS.hg38.vardict.100bp_padded.bed")

########################## Calculate and Set Up Data ###########################
i=-1
for dataset in ${datasets[@]}
  do
  let i=${i}+1
  echo "Processing dataset: ${dataset}"
  Rscript analyses/snv-callers/scripts/01-calculate_vaf_tmb.R \
  -m data/pbta-snv-${dataset}.vep.maf.gz \
  -b data/${wgs_files[$i]} \
  -l ${dataset}
  done
