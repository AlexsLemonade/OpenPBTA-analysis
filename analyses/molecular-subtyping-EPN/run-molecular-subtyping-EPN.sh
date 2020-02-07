#!/bin/bash


# Run `00-subsetting-files-for-EPN.py`  to subset gene expression data and 
# Run `01-make_notebook_RNAandDNA.py` to match RNA and DNA Biospecimen ID's based on corresponding matching  sample_ID from pbta-histologies.tsv  file and 
# Run `02_ependymoma_generate_all_data.py` to combine data  from various output files for EPN samples that can be used to categorize samples under different EPN subtypes 

set -e 
set -o pipefail 

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# This option controls whether on not the step that generates the EPN only
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

# Define needed files
HISTOLOGIES=../../data/pbta-histologies.tsv
FULL_EXPRESSION=../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
SUBSET_EXPRESSION=epn-subset/epn-pbta-gene-expression-rsem-fpkm-collapsed.stranded.tsv.gz
NOTEBOOK=results/EPN_molecular_subtype.tsv
GISTIC=../../data/pbta-cnv-consensus-gistic.zip
GISTIC_FILE=pbta-cnv-consensus-gistic/broad_values_by_arm.txt
GSVA=../gene-set-enrichment-analysis/results/gsva_scores_stranded.tsv
FUSION=../fusion-summary/results/fusion_summary_ependymoma_foi.tsv
BREAKPOINTS=../chromosomal-instability/breakpoint-data/union_of_breaks_densities.tsv

OUTPUT=results/EPN_all_data.tsv

# make the subset and results directory if they don't exist
mkdir -p epn-subset
mkdir -p results

if [ "$SUBSET" -gt "0" ]; then
  echo "Subsetting  for CI"
  Rscript 00-subset-for-EPN.R -i $HISTOLOGIES -e $FULL_EXPRESSION -o $SUBSET_EXPRESSION 
fi

echo "Generating analyses/molecular-subtyping-EPN/results/EPN_molecular_subtype.tsv that maps DNA and RNA ID's"
python3 01-make_notebook_RNAandDNA.py -i $HISTOLOGIES -o $NOTEBOOK

echo  "Generating analyses/molecular-subtyping-EPN/results/EPN_all_data.tsv  that has all the relevant data needed for subtyping"
python3 02_ependymoma_generate_all_data.py \
    --notebook $NOTEBOOK \
    --gistic $GISTIC \
    --subfile-gistic $GISTIC_FILE \
    --gsva $GSVA \
    --expression $SUBSET_EXPRESSION \
    --fusion $FUSION \
    --breakpoints $BREAKPOINTS \
    --outfile $OUTPUT

