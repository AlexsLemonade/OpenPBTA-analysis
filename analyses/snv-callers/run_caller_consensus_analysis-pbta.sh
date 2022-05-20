#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: Run a consensus analysis for PBTA of SNV callers

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# The sqlite database made from the callers will be called:
dbfile=../../scratch/snv-callers/snv_db.sqlite

# Designate output file
consensus_file=results/consensus/pbta-snv-consensus-mutation.maf.tsv

# BED and GTF file paths
cds_file=../../scratch/gencode.v27.primary_assembly.annotation.bed
wgs_bed=../../scratch/intersect_strelka_mutect_WGS.bed

# Set a default for the VAF filter if none is specified
vaf_cutoff=${OPENPBTA_VAF_CUTOFF:-0}

# If running during data generation, we want to use the BASE histologies file
run_for_release=${OPENPBTA_BASE_RELEASE:-0}
if [[ "$run_for_release" -eq "0" ]]
then
   histologies_file=../../data/pbta-histologies.tsv
else
   histologies_file=../../data/pbta-histologies-base.tsv
fi

# Unless told to run the plots, the default is to skip them
# To run plots, set OPENPBTA_PLOTS to 1 or more
run_plots_nb=${OPENPBTA_PLOTS:-0}

################################ Set Up Database ################################
echo "Setting up Database"
python3 scripts/01-setup_db.py \
  --db-file $dbfile \
  --strelka-file ../../data/pbta-snv-strelka2.vep.maf.gz \
  --mutect-file ../../data/pbta-snv-mutect2.vep.maf.gz \
  --lancet-file ../../data/pbta-snv-lancet.vep.maf.gz \
  --vardict-file ../../data/pbta-snv-vardict.vep.maf.gz \
  --meta-file $histologies_file

##################### Merge callers' files into total files ####################
echo "Merging callers"
Rscript scripts/02-merge_callers.R \
  --db_file $dbfile \
  --output_file $consensus_file \
  --vaf_filter $vaf_cutoff \
  --overwrite

########################## Add consensus to db ################################
echo "Adding consensus to database"
python3 scripts/01-setup_db.py \
  --db-file $dbfile \
  --consensus-file $consensus_file

############# Create intersection BED files for TMB calculations ###############
# Make All mutations BED files
echo "Making intersection bed files"
bedtools intersect \
  -a ../../data/WGS.hg38.strelka2.unpadded.bed \
  -b ../../data/WGS.hg38.mutect2.vardict.unpadded.bed \
  > $wgs_bed

#################### Make coding regions file
# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
echo "Making CDS bed file"
gunzip -c ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge  \
  > $cds_file

######################### Calculate consensus TMB ##############################
echo "Calculating TMB"
Rscript scripts/03-calculate_tmb.R \
  --db_file $dbfile \
  --output results/consensus \
  --metadata $histologies_file \
  --coding_regions $cds_file \
  --overwrite \
  --nonsynfilter_maf

########################## Compress consensus file #############################

gzip $consensus_file

############################# Comparison Plots #################################
if [ "$run_plots_nb" -gt "0" ]
then
 echo "Making comparison plots"
 Rscript -e "rmarkdown::render('compare_snv_callers_plots.Rmd', clean = TRUE)"
fi
