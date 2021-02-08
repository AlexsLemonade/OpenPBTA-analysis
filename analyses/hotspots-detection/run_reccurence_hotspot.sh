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


# The sqlite database made from the callers will be called:
dbfile=../../scratch/snv_db.sqlite

################################ Set Up Database ################################
##### Updated form snv-caller module to include all columns for mutect and vardict

#python3 00-setup_db.py \
#  --db-file $dbfile \
#  --strelka-file ../../data/pbta-snv-strelka2.vep.maf.gz \
#  --mutect-file ../../data/pbta-snv-mutect2.vep.maf.gz \
#  --lancet-file ../../data/pbta-snv-lancet.vep.maf.gz \
#  --vardict-file ../../data/pbta-snv-vardict.vep.maf.gz \


# find reccurence
Rscript -e "rmarkdown::render('01-reccurence-hotspot-overlap.Rmd')"
