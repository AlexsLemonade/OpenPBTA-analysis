# C. Bethell and C. Savonen for CCDL 2019
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

XYFLAG=${OPENPBTA_XY:-1}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Prep the CNVkit data
Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Run annotation step for CNVkit
Rscript --vanilla 03-prepare-cn-file.R \
  --cnv_file ../../scratch/cnvkit_with_status.tsv \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "cnvkit_annotated_cn" \
  --seg

# Run annotation step for ControlFreeC
Rscript --vanilla 03-prepare-cn-file.R \
  --cnv_file ../../data/pbta-cnv-controlfreec.tsv.gz \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "controlfreec_annotated_cn" \
  --xy $XYFLAG \
  --controlfreec

# Compare to expression data
# Rscript --vanilla 04-rna-expression-validation.R
