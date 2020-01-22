#!/bin/bash
# CCDL for ALSF 2020
# Candace L Savonen
# 
###################### Create intersection BED files ###########################
# WGS effectively surveyed BED file
# Note strelka2's BED file is being used here because all it includes the whole genome
# Then we are subtracting the ranges from the unsurveyed regions which were identified in
# 
surveyed_wgs=data/WGS.hg38.strelka2.unpadded.bed

# TODO: when the unsurveyed regions of the genome BED file is prepared, replace this
# <unsurveyed_regions.bed> with the appropriate file and uncomment this section. 
# Also change the above file name to a new file name. 

# bedtools subtract \
#  -a data/WGS.hg38.strelka2.unpadded.bed \
#  -b <unsurveyed_regions.bed> > $surveyed_wgs
  
# WXS effectively surveyed BED file
# TODO: potentially we should be using an unpadded BED file here to get a more accurate estimate of the surveyed genome. 
# TODO: we may also want to subtract <unsurveyed_regions.bed> regions from this BED regions. 
surveyed_wxs=data/WXS.hg38.100bp_padded.bed  

############################ Run setup data script #############################
# TODO: change CNVKit file to consensus CNV file
Rscript analyses/chromosomal-instability/00-setup-breakpoint-data.R \
  --cnv_seg data/pbta-cnv-cnvkit.seg.gz \
  --sv data/pbta-sv-manta.tsv.gz \
  --metadata data/pbta-histologies.tsv \
  --output analyses/chromosomal-instability/breakpoint-data \
  --surveyed_wgs $surveyed_wgs \
  --surveyed_wxs $surveyed_wxs
  
######################### Chromosomal Instability Plots ########################
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/01-plot-chromosomal-instability.Rmd', 
                              clean = TRUE)"
