# Chante Bethell for CCDL 2019
# Run 01-prepare-cn-file.R
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
Rscript --vanilla -e "rmarkdown::render('00-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Run annotation step for CNVkit
Rscript --vanilla 01-prepare-cn-file.R \
  --cnv_file ../../scratch/cnvkit_with_status.tsv \
  --gtf_file ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "cnvkit_annotated_cn" \
  --cnvkit

# Run annotation step for ControlFreeC
Rscript --vanilla 01-prepare-cn-file.R \
  --cnv_file ../../data/pbta-cnv-controlfreec.tsv.gz \
  --gtf_file ../../data/gencode.v27.primary_assembly.annotation.gtf.gz  \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "controlfreec_annotated_cn" \
  --xy $XYFLAG \
  --controlfreec

# gzip the four files in the results folder, overwriting without prompt
gzip -f results/cnvkit_annotated_cn_autosomes.tsv
gzip -f results/controlfreec_annotated_cn_autosomes.tsv
gzip -f results/cnvkit_annotated_cn_x_and_y.tsv

if [ "$XYFLAG" -gt "0" ]; then
  gzip -f results/controlfreec_annotated_cn_x_and_y.tsv
fi

# Compare to expression data
Rscript --vanilla 02-rna-expression-validation.R
