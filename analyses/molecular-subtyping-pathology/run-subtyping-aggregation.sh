#!/bin/bash
#
# J. Taroni for ALSF CCDL 2020
#
# Run aggregation of molecular subtyping/reclassification results and the
# incorporation of pathology feedback

set -e
set -o pipefail

# We're using this to tie the clinical file to a specific release when this
# is not run in CI
IS_CI=${OPENPBTA_TESTING:-0}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Run the first notebook that compiles all the results from other modules into
# a single table
Rscript -e "rmarkdown::render('01-compile-subtyping-results.Rmd', params=list(is_ci = ${IS_CI}), clean = TRUE)"

# Recoding ACP samples
Rscript -e "rmarkdown::render('pathology-subtyping-craniopharyngioma.Rmd', clean = TRUE)"

# Run the second notebook to incorporate clinical review to the compiled subtyping
Rscript -e "rmarkdown::render('02-incorporate-clinical-feedback.Rmd', clean = TRUE)"

# Run the third notebook that incorporates pathology feedback into final labels
Rscript -e "rmarkdown::render('03-incorporate-pathology-feedback.Rmd', params=list(is_ci = ${IS_CI}), clean = TRUE)"

# Run the meningioma pathology-free-text based subtyping step
Rscript -e "rmarkdown::render('pathology_free_text-subtyping-meningioma.Rmd', clean = TRUE)"

# Glialneuronal tumors
Rscript -e "rmarkdown::render('pathology-harmonized-diagnosis-glialneuronal-tumors.Rmd', clean = TRUE)"

# Choroid plexus papilloma
Rscript -e "rmarkdown::render('pathology-subtyping-choroid-plexus-papilloma.Rmd', clean = TRUE)"

