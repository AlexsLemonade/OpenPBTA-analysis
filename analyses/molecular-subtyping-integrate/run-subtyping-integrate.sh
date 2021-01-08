#!/bin/bash
#
# Krutika Gaonkar D3b
#
# Integrate molecular subtyping results

set -e
set -o pipefail

# We're using this to tie the clinical file to a specific release when this
# is not run in CI
IS_CI=${OPENPBTA_TESTING:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Add molecular_subtype,	integrated_diagnosis, short_histology,	
# broad_histology and	Notes from `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`.
# Check for changes and duplicates and save a final file
Rscript -e "rmarkdown::render('01-integrate-subtyping.Rmd',params=list(is_ci = ${IS_CI}),clean = TRUE)"
