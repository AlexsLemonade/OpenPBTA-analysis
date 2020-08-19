#!/bin/bash
# Module author: Komal S. Rathi
# 2020

# This script runs the steps for molecular subtyping of Medulloblastoma samples

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# create results directory if it doesn't already exist
mkdir -p results
mkdir -p plots

# run MB molecular subtyping (no batch correction)
Rscript --vanilla 01-classify-mb.R \
--polyaexprs ../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
--strandedexprs ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
--clin ../../data/pbta-histologies.tsv \
--batch_col NA \
--hist_column integrated_diagnosis \
--method medullo-classifier \
--outputfile results/mb-molecular-subtypes-medullo-classifier.rds

# merge expected and observed outputs
Rscript 02-compare-classes.R \
--observed_class results/mb-molecular-subtypes-medullo-classifier.rds \
--expected_class input/expected_class.rds \
--outputfile results/comparison-medullo-classifier.rds