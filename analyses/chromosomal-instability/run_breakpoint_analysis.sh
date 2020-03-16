#!/bin/bash
# CCDL for ALSF 2020
# Candace L Savonen
#
# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# Need to adjust minimum samples to plot if in CI:
IS_CI=${OPENPBTA_TESTING:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Data and scratch directories
data_dir=../../data
scratch_dir=../../scratch

# To avoid conflicts, we're going to create a folder within scratch specifically
# for this module
module_scratch_dir=${scratch_dir}/chromosomal-instability
mkdir -p ${module_scratch_dir}

# We use information from the copy_number_consensus_call module such as the
# uncalled regions
cn_consensus_dir=../copy_number_consensus_call

###################### Create intersection BED files ###########################
# WGS effectively surveyed BED file
surveyed_wgs=${module_scratch_dir}/cnv_surveyed_wgs.bed

# TODO: when the unsurveyed regions of the genome BED file is add to the data
# release, update this file path
# Note strelka2's BED file is being used here because all it includes the whole genome
 bedtools subtract \
  -a ${data_dir}/WGS.hg38.strelka2.unpadded.bed \
  -b ${cn_consensus_dir}/ref/cnv_excluded_regions.bed > $surveyed_wgs

# WXS effectively surveyed BED file
# Currently there is no WXS CNV data in the data release, so this WXS BED file is 
# just a place holder for if that happens and is not actually used. 
# TODO: If WXS data is added, we may also want to subtract <unsurveyed_regions.bed> 
# regions from this BED regions as well (like is done above with WGS)
surveyed_wxs=${data_dir}/WXS.hg38.100bp_padded.bed

############################ Run setup data script #############################
Rscript 00-setup-breakpoint-data.R \
  --cnv_seg ${data_dir}/pbta-cnv-consensus.seg.gz \
  --sv ${data_dir}/pbta-sv-manta.tsv.gz \
  --metadata ${data_dir}/pbta-histologies.tsv \
  --uncalled ${cn_consensus_dir}/results/uncalled_samples.tsv \
  --output breakpoint-data \
  --surveyed_wgs $surveyed_wgs \
  --surveyed_wxs $surveyed_wxs \
  --gap 5 \
  --drop_sex

######################### Localization calculations ############################
Rscript -e "rmarkdown::render('01-localization-of-breakpoints.Rmd',
                              clean = TRUE)"

######################### Chromosomal Instability Plots ########################
# Circos plots examples:
Rscript -e "rmarkdown::render('01b-visualization-cnv-sv.Rmd',
                              clean = TRUE)"
# Heatmaps:
Rscript -e "rmarkdown::render('02a-plot-chr-instability-heatmaps.Rmd',
                              clean = TRUE)"
# Histology plots:
if [ $IS_CI -gt 0 ]
then
  MIN_SAMPLES=0
else
  MIN_SAMPLES=5
fi
echo $MIN_SAMPLES
Rscript -e "rmarkdown::render('02b-plot-chr-instability-by-histology.Rmd',
                              clean = TRUE, params = list(min_samples=${MIN_SAMPLES}))"
