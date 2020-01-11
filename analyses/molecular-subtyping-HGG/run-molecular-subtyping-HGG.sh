#!/bin/bash

# Chante Bethell for CCDL 2020
#
# Run `01-HGG-molecular-subtyping-defining-lesions.Rmd` and
# `02-HGG-molecular-subtyping-subset-files.R` if needed.

set -e
set -o pipefail

# This option controls whether on not the step that generates the HGG only
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

# cds gencode bed file is used by other analyses where mutation data is
# filtered to only coding regions
exon_file="../../scratch/gencode.v27.primary_assembly.annotation.bed"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run the first script in this module that reclassifies high-grade gliomas
Rscript -e "rmarkdown::render('01-HGG-molecular-subtyping-defining-lesions.Rmd', clean = TRUE)"

# Run the second script in this module that subset files using the samples in the output
# file generated with `01-HGG-molecular-subtyping-defining-lesions.Rmd`.
if [ "$SUBSET" -gt "0" ]; then
  Rscript --vanilla 02-HGG-molecular-subtyping-subset-files.R
fi

#### Mutation data -------------------------------------------------------------

# if the cds gencode bed file is not available from another analysis, generate
# it here
if [ ! -f "$exon_file" ]; then
  gunzip -c ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
    | awk '$3 ~ /CDS/' \
    | convert2bed --do-not-sort --input=gtf - \
    > $exon_file
if

