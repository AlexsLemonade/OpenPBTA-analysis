# C. Bethell and C. Savonen for CCDL 2019
# Run 00-add-ploidy-cnvkit.Rmd, 01-prepare-cn-file.R, and 02-rna-expression-validation.R
# sequentially.
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# This cds file is filtered below to contain coding sequences only, converted
# to a bed file, and saved in the project's `scratch` directory.
XYFLAG=${OPENPBTA_XY:-1}
cds_file=../../scratch/gencode.v27.primary_assembly.annotation.bed

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up cds file
gunzip -c ../../data/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > $cds_file

# Prep the CNVkit data
Rscript --vanilla -e "rmarkdown::render('00-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Run annotation step for CNVkit
Rscript --vanilla 01-prepare-cn-file.R \
  --cnv_file ../../scratch/cnvkit_with_status.tsv \
  --cds_file $cds_file \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "cnvkit_annotated_cn" \
  --cnvkit \
  --gistic

# Run annotation step for ControlFreeC
Rscript --vanilla 01-prepare-cn-file.R \
  --cnv_file ../../data/pbta-cnv-controlfreec.tsv.gz \
  --cds_file $cds_file \
  --metadata ../../data/pbta-histologies.tsv \
  --filename_lead "controlfreec_annotated_cn" \
  --xy $XYFLAG \
  --controlfreec \
  --gistic

# Compare to expression data
Rscript --vanilla 02-rna-expression-validation.R
