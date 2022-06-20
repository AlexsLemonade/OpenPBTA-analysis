#!/bin/bash
# Module author: Komal S. Rathi
# 2019

# This script runs the steps for immune deconvolution in PBTA histologies using xCell. 
# xCell is the most comprehensive deconvolution method (with the largest number of cell types) and widely used in literature vs other deconvolution methods.
# Reference for benchmarking between xCell and other methods: PMID: 31641033

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

# Shared input files
POLYA='../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds'
STRANDED='../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds'
CLIN='../../data/pbta-histologies.tsv'

# In CI we'll run an abbreviated version of the figures script due to insufficient testing data
# needed for making some of the plots
ABBREVIATED_IMMUNE=${OPENPBTA_QUICK_IMMUNE:-0}

# generate deconvolution output for poly-A and stranded datasets using xCell
Rscript --vanilla 01-immune-deconv.R \
--polyaexprs $POLYA \
--strandedexprs $STRANDED \
--method 'xcell' \
--outputfile 'results/xcell_deconv-output.rds' 

# generate deconvolution output for poly-A and stranded datasets using quanTIseq
Rscript --vanilla 01-immune-deconv.R \
--polyaexprs $POLYA \
--strandedexprs $STRANDED \
--method 'quantiseq' \
--outputfile 'results/quantiseq_deconv-output.rds' 


echo "Deconvolution finished."


# Perform visualization
Rscript -e "rmarkdown::render('02-visualize_quantiseq.Rmd', clean = TRUE, params = list(is_ci = ${ABBREVIATED_IMMUNE}))"

