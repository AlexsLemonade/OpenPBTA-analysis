#!/bin/bash
# CCDL for ALSF 2020
# Candace L Savonen
# 
###################### Create intersection BED files ###########################
# WGS effectively surveyed BED file
surveyed_wgs=scratch/cnv_surveyed_wgs.bed

# TODO: when the unsurveyed regions of the genome BED file is add to the data 
# release, update this file path
# Note strelka2's BED file is being used here because all it includes the whole genome
 bedtools subtract \
  -a data/WGS.hg38.strelka2.unpadded.bed \
  -b analyses/copy_number_consensus_call/ref/cnv_excluded_regions.bed > $surveyed_wgs
  
# WXS effectively surveyed BED file
# TODO: potentially we should be using an unpadded BED file here to get a more accurate estimate of the surveyed genome. 
# TODO: we may also want to subtract <unsurveyed_regions.bed> regions from this BED regions. 
surveyed_wxs=data/WXS.hg38.100bp_padded.bed  

############################ Run setup data script #############################
# TODO: update CNV consensus file path with v14
Rscript analyses/chromosomal-instability/00-setup-breakpoint-data.R \
  --cnv_seg analyses/copy_number_consensus_call/results/pbta-cnv-consensus.seg.gz \
  --sv data/pbta-sv-manta.tsv.gz \
  --metadata data/pbta-histologies.tsv \
  --output analyses/chromosomal-instability/breakpoint-data \
  --surveyed_wgs $surveyed_wgs \
  --surveyed_wxs $surveyed_wxs \
  --gap 5 \
  --drop_sex
  
######################### Localization calculations ############################
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/01-localization-of-breakpoints.Rmd', 
                              clean = TRUE)"
                              
######################### Chromosomal Instability Plots ########################   
# Circos plots examples:
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/01b-visualization-cnv-sv.Rmd', 
                              clean = TRUE)"
# Heatmaps:
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/02a-plot-chr-instability-heatmaps.Rmd', 
                              clean = TRUE)"
