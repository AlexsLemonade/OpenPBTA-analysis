#!/bin/bash

set -e
set -o pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# This option controls whether on not the step that generates the LGAT only mutation 
# files gets run -- it will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}


if [ "$SUBSET" -gt "0" ]; then
# subset bu SNV 
Rscript -e "rmarkdown::render('01-subset-files-for-LGAT.Rmd')"
# subset by Fusion
Rscript -e "rmarkdown::render('02-subset-fusion-files-LGAT.Rmd')"
# subset by CNV
Rscript -e "rmarkdown::render('03-subset-cnv-files-LGAT.Rmd')"
fi


# compile subtypes
Rscript -e "rmarkdown::render('04-LGAT-compile-subtypes.Rmd')"
