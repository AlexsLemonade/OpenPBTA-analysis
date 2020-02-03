# C. Bethell and C. Savonen for CCDL 2019
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

XYFLAG=${OPENPBTA_XY:-1}
# currently, we use the consensus SEG file committed to the repository
# this means it is *not* subset and therefore requires too much RAM
# to run the annotation step
RUNCONSENSUS=${OPENPBTA_CONSENSUS:-1}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
histologies_file=${data_dir}/pbta-histologies.tsv
gtf_file=${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz

# Prep the CNVkit data
Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

# Run annotation step for CNVkit
Rscript --vanilla 03-prepare-cn-file.R \
  --cnv_file ${scratch_dir}/cnvkit_with_status.tsv \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "cnvkit_annotated_cn" \
  --seg

# Run annotation step for ControlFreeC
Rscript --vanilla 03-prepare-cn-file.R \
  --cnv_file ${data_dir}/pbta-cnv-controlfreec.tsv.gz \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "controlfreec_annotated_cn" \
  --xy $XYFLAG \
  --controlfreec


if [ "$RUNCONSENSUS" -gt "0"]; then
# Run annotation step for consensus file
  Rscript --vanilla 03-prepare-cn-file.R \
    --cnv_file ${scratch_dir}/consensus_seg_with_status.tsv \
    --gtf_file $gtf_file \
    --metadata $histologies_file \
    --filename_lead "consensus_seg_annotated_cn" \
    --seg \
    --xy $XYFLAG
fi

# Loop over all annotated cn files in the results directory for each of the
# collapsed expression files (polyA and stranded)
FILES=results/*
for f in $FILES
do
  echo "Plotting $f ..."
  # We want to extract the `filename_lead` from the name of the file in the
  # results directory
  filename_lead=${f%%.*}
  filename_lead=${filename_lead#*/}
  filename_lead=${filename_lead}_polya
  Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file $f \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
  --independent_specimens_file ../../data/independent-specimens.wgswxs.primary.tsv \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead $filename_lead

  filename_lead=${f%%.*}
  filename_lead=${filename_lead#*/}
  filename_lead=${filename_lead}_stranded
  Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file $f \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
  --independent_specimens_file ../../data/independent-specimens.wgswxs.primary.tsv \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead $filename_lead
done

