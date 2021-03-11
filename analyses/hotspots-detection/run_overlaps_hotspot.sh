#!/bin/bash

# K S Gaonkar
set -e
set -o pipefail


# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

cancer_hotspot_file=input/hotspots_v2.xls
genomic_site_hotspot_file=input/tert_promoter_hotspots.tsv


# The sqlite database made from the callers will be called:
dbfile=../../scratch/snv_db.sqlite

################################ Set Up Database ################################
##### this script from snv-caller module is being used to set up a database with tables for each caller

#python3 ../snv-callers/scripts/01-setup_db.py \
#  --db-file $dbfile --overwrite \
#  --strelka_file ../../data/pbta-snv-strelka2.vep.maf.gz
#  --lancet_file ../../data/pbta-snv-lancet.vep.maf.gz
#  --mutect-file ../../data/pbta-snv-mutect2.vep.maf.gz \
#  --vardict-file ../../data/pbta-snv-vardict.vep.maf.gz \
################################################################################


# find hotspot overlaps in strelka2
#Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-strelka2.vep.maf.gz \
#        --caller strelka2 \
#        --cancer_hotspot_file $cancer_hotspot_file \
#        --genomic_site_hotspot_file $genomic_site_hotspot_file 

# find hotspot overlaps in mutect2
#Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-mutect2.vep.maf.gz \
#        --caller mutect2 \
#        --cancer_hotspot_file $cancer_hotspot_file \
#        --genomic_site_hotspot_file $genomic_site_hotspot_file

# find hotspot overlaps in lancet
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-lancet.vep.maf.gz \
        --caller lancet \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file


# find hotspot overlaps in vardict
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-vardict.vep.maf.gz \
        --caller vardict \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file

#Rscript -e "rmarkdown::render('01-combine-maf.Rmd')"
