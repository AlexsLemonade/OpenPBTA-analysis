#!/bin/bash
set -eo pipefail

# set the analyses folder as workdir variable
cd "$(dirname "${BASH_SOURCE[0]}")"

# retrieve the TCGA's exome capture kit from GDC's file API endpoint
# this will generate tcga-capture_kit-info.tsv file under the result folder
python3 scripts/get-tcga-capture_kit.py

# download all BED from tcga-capture_kit-info.tsv and add chr prefix
bash scripts/prepare-tcga-capture_kit.sh

# Convert BED from hg19 to Gh38
bash scripts/CrossMap.sh

## get somatic mutation counts for intersection regions
bash scripts/intersect-bed-maf.sh

## boxplot from ggplot and saved to plots folder
Rscript scripts/boxplot.R ../../scratch/somatic-count_with-histologies.tsv plots/

