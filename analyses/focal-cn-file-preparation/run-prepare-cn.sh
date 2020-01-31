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
  
# Compare annotated CNVkit autosome output to polyA expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/cnvkit_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.polya.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "cnvkit_annotated_cn_autosomes_polya"

# Compare annotated CNVkit autosome output to stranded expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/cnvkit_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "cnvkit_annotated_cn_autosomes_stranded"
  
# Compare annotated Controlfreec autosome output to polyA expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/controlfreec_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.polya.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "controlfreec_annotated_cn_autosomes_polya"

# Compare annotated Controlfreec autosome output to stranded expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/controlfreec_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "controlfreec_annotated_cn_autosomes_stranded"

if [ "$RUNCONSENSUS" -gt "0"]; then
# Compare annotated consensus autosome output to polyA expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/consensus_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.polya.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "consensus_annotated_cn_autosomes_polya"

# Compare annotated consensus autosome output to stranded expression data 
Rscript --vanilla rna-expression-validation.R \
  --annotated_cnv_file results/consensus_annotated_cn_autosomes.tsv.gz \
  --expression_file ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata $histologies_file \
  --goi_list ../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
  --filename_lead "consensus_annotated_cn_autosomes_stranded"

fi
