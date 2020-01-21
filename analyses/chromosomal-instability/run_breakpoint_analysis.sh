#!/bin/bash
# CCDL for ALSF 2020
# Candace L Savonen
# 
###################### Create intersection BED files ###########################
# WGS effectively surveyed BED file
# Note strelka2's BED file is being used here because all it includes the whole genome
# Then we are subtracting the ranges from the unsurveyed regions which were identified in ...
bedtools subtract \
  -a data/WGS.hg38.strelka2.unpadded.bed \
  -b data/unsurveyed_regions.bed > WGS_effectively_surveyed.bed
  
# WXS effectively surveyed BED file
data/WXS.hg38.100bp_padded.bed  

############################ Run setup data script #############################
Rscript analyses/chromosomal-instability/00-setup-breakpoint-data.R \
  --cnv_seg data/pbta-cnv-cnvkit.seg.gz \
  --sv data/pbta-sv-manta.tsv.gz \
  --metadata data/pbta-histologies.tsv \
  --surveyed_wgs  data/WGS.hg38.strelka2.unpadded.bed \
  --surveyed_wxs data/WXS.hg38.100bp_padded.bed
  
######################### Chromosomal Instability Plots ########################
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/plot-chromosomal-instability.Rmd', 
                              clean = TRUE)"
                              
############################ Recurrent Genes Analysis ##########################
# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c data/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > scratch/exon_gencode.v27.primary_assembly.annotation.bed

# Run the recurrent gene notebook 
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/find-recurrent-breakpoint-genes.Rmd', 
                              clean = TRUE)"
