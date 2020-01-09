#!/bin/bash
# Module author: Komal S. Rathi
# 2019

# This script runs the steps for immune deconvolution in PBTA histologies
# It uses xCell and one of the three methods:
# mcp_counter, cibersort_abs, timer

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
# generate deconvolution output for poly-A and stranded datasets using xCell and the second specified method
# we will use CIBERSORT for the paper but because the dependent scripts are not publicly accesible, we will use the next best method i.e. MCP-counter as default for testing purpose
# Reference for benchmarking between xCell, CIBERSORT and MCP-counter: PMID: 31641033

DECONV_METHOD=${OPENPBTA_DECONV_METHOD:-"cibersort_abs"}
echo "Deconv method: $DECONV_METHOD"
if [ "$DECONV_METHOD" == "cibersort_abs" ]
then
	CIBERSORT_BIN="CIBERSORT.R"
	CIBERSORT_MAT="LM22.txt"
else
	CIBERSORT_BIN="NA"
	CIBERSORT_MAT="NA"
fi
echo "$CIBERSORT_BIN"
echo "$CIBERSORT_MAT"

Rscript --vanilla 01-immune-deconv.R \
--polyaexprs ../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
--strandedexprs ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
--clin ../../data/pbta-histologies.tsv \
--method $DECONV_METHOD \
--cibersortbin $CIBERSORT_BIN \
--cibersortgenemat $CIBERSORT_MAT \
--outputfile results/deconv-output.RData

echo "Deconvolution finished..."
echo "Create summary plots"

# Now, run the script to generate correlation plots between xCell and the second method 
# Also generates corresponding heatmaps of average normalized immune scores
Rscript --vanilla 02-summary-plots.R \
--input results/deconv-output.RData \
--output plots

