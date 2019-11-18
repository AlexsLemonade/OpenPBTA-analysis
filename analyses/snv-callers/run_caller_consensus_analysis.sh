#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose:Run an consensus analysis of SNV callers

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# The sqlite database made from the callers will be called:
dbfile=scratch/testing_snv_db.sqlite

# Designate output file 
consensus_file=analyses/snv-callers/results/consensus/consensus_snv.maf.tsv

# Reference file paths
cosmic=analyses/snv-callers/ref_files/brain_cosmic_variants_coordinates.tsv
annot_rds=analyses/snv-callers/ref_files/hg38_genomic_region_annotation.rds

# Set a default for the VAF filter if none is specified
vaf_cutoff=${OPENPBTA_VAF_CUTOFF:-0}

############################ Set Up Reference Files ############################
# The original COSMIC file is obtained from: https://cancer.sanger.ac.uk/cosmic/download
# These data are available if you register. The full, unfiltered somatic mutations
# file CosmicMutantExport.tsv.gz for grch38 is used here.
Rscript analyses/snv-callers/scripts/00-set_up.R \
  --annot_rds $annot_rds \
  --cosmic_og scratch/CosmicMutantExport.tsv.gz \
  --cosmic_clean $cosmic

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
  
######################### Calculate consensus TMB ##############################
Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
  --consensus analyses/snv-callers/results/consensus/consensus_snv.maf.tsv \
  --output analyses/snv-callers/results/consensus \
  --metadata data/pbta-histologies.tsv \
  --bed_wgs data/WGS.hg38.strelka2.unpadded.bed \
  --bed_wxs data/WXS.hg38.100bp_padded.bed \
  --overwrite
