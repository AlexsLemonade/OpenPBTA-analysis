#!/bin/bash


# Run `00-subsetting-files-for-EPN.py`  to subset gene expression data and 
# Run `01-make_notebook_RNAandDNA.py` to match RNA and DNA Biospecimen ID's based on corresponding matching  sample_ID from pbta-histologies.tsv  file and 
# Run `02_ependymoma_generate_all_data.py` to combine data  from various output files for EPN samples that can be used to categorize samples under different EPN subtypes 

set -e 
set -o pipefail 


# This option controls whether on not the step that generates the EPN only
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

if [ "$SUBSET" -gt "0" ]; then
  echo "Subsettibg  for CI"
  python3 analyses/molecular-subtyping-EPN/00-subsetting-files-for-EPN.py
fi

echo "Generating analyses/molecular-subtyping-EPN/results/EPN_molecular_subtype.tsv that maps DNA and RNA ID's"
python3 analyses/molecular-subtyping-EPN/01-make_notebook_RNAandDNA.py

echo  "Generating analyses/molecular-subtyping-EPN/results/EPN_all_data.tsv  that has all the relevant data needed for subtyping"
python3 analyses/molecular-subtyping-EPN/02_ependymoma_generate_all_data.py

