#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: Run a consensus analysis for PBTA of SNV callers

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# The sqlite database made from the callers will be called:
dbfile=scratch/snv_db.sqlite

# Designate output file 
consensus_file=analyses/snv-callers/results/consensus/pbta-snv-consensus-mutation.maf.tsv

# BED and GTF file paths
cds_file=scratch/gencode.v27.primary_assembly.annotation.bed
all_mut_wgs_bed=scratch/intersect_strelka_mutect_WGS.bed
all_mut_wxs_bed=data/WXS.hg38.100bp_padded.bed
coding_wgs_bed=scratch/intersect_cds_strelka_mutect_WGS.bed
coding_wxs_bed=scratch/intersect_cds_strelka_mutect_WXS.bed

# Set a default for the VAF filter if none is specified
vaf_cutoff=${OPENPBTA_VAF_CUTOFF:-0}

# Unless told to run the plots, the default is to skip them
# To run plots, set OPENPBTA_PLOTS to 1 or more
run_plots_nb=${OPENPBTA_PLOTS:-0}

################################ Set Up Database ################################
python3 analyses/snv-callers/scripts/01-setup_db.py \
  --db-file $dbfile \
  --strelka-file data/pbta-snv-strelka2.vep.maf.gz \
  --mutect-file data/pbta-snv-mutect2.vep.maf.gz \
  --lancet-file data/pbta-snv-lancet.vep.maf.gz \
  --vardict-file data/pbta-snv-vardict.vep.maf.gz \
  --meta-file data/pbta-histologies.tsv

##################### Merge callers' files into total files ####################
Rscript analyses/snv-callers/scripts/02-merge_callers.R \
  --db_file $dbfile \
  --output_file $consensus_file \
  --vaf_filter $vaf_cutoff \
  --overwrite

########################## Add consensus to db ################################
python3 analyses/snv-callers/scripts/01-setup_db.py \
  --db-file $dbfile \
  --consensus-file $consensus_file

############# Create intersection BED files for TMB calculations ###############
# Make All mutations BED files
bedtools intersect \
  -a data/WGS.hg38.strelka2.unpadded.bed \
  -b data/WGS.hg38.mutect2.vardict.unpadded.bed \
  > $all_mut_wgs_bed

#################### Make coding regions file 
# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c data/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge  \
  > $cds_file
  
##################### Make WGS coding BED file  
# Make WGS coding BED file for strelka
bedtools intersect \
  -a data/WGS.hg38.strelka2.unpadded.bed \
  -b $cds_file \
  > scratch/wgs_coding_strelka.bed

# Make WGS coding BED file for mutect
bedtools intersect \
  -a data/WGS.hg38.mutect2.vardict.unpadded.bed \
  -b $cds_file \
  > scratch/wgs_coding_mutect.bed

# Intersect the mutect and strelka coding beds into one
bedtools intersect \
  -a scratch/wgs_coding_mutect.bed \
  -b scratch/wgs_coding_strelka.bed \
   > scratch/wgs_coding_strelka_mutect.bed

# Merge these ranges into one
sort -k 1,1 -k 2,2n scratch/wgs_coding_strelka_mutect.bed | 
bedtools merge -i stdin > $coding_wgs_bed
   
##################### Make WXS coding BED file
# Intersect coding and WXS ranges, sort and merge 
bedtools intersect \
  -a data/WXS.hg38.100bp_padded.bed  \
  -b $cds_file |
sort -k 1,1 -k 2,2n -i stdin |
bedtools merge -i stdin > $coding_wxs_bed

######################### Calculate consensus TMB ##############################
Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
  --db_file $dbfile \
  --output analyses/snv-callers/results/consensus \
  --metadata data/pbta-histologies.tsv \
  --all_bed_wgs $all_mut_wgs_bed \
  --all_bed_wxs $all_mut_wxs_bed \
  --coding_bed_wgs $coding_wgs_bed \
  --coding_bed_wxs $coding_wxs_bed \
  --overwrite
 
########################## Compress consensus file #############################

gzip $consensus_file

############################# Comparison Plots #################################
if [ "$run_plots_nb" -gt "0" ]
then
 Rscript -e "rmarkdown::render('analyses/snv-callers/compare_snv_callers_plots.Rmd', clean = TRUE)"
fi
