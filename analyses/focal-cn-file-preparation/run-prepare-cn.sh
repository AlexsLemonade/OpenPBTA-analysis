# C. Bethell and C. Savonen for CCDL 2019
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# Run original files - will not by default
RUN_ORIGINAL=${RUN_ORIGINAL:-0}

# Run testing files for circle CI - will not by default
IS_CI=${OPENPBTA_TESTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
results_dir=../../analyses/focal-cn-file-preparation/results
histologies_file=${data_dir}/pbta-histologies.tsv
gtf_file=${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz
goi_file=../../analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt
independent_specimens_file=${data_dir}/independent-specimens.wgswxs.primary.tsv

# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

# Run snakemake script implementing `bedtools coverage` for each sample bed file in
# `scratch/cytoband_status` -- these files are generated in
# `02-add-ploidy-consensus.Rmd`
# currently runs 10 jobs in parallel, which should be fine for most implementations
snakemake -j 10 --snakefile run-bedtools.snakemake

# Determine the dominant status for each chromosome arm and compare with GISTIC's arm status calls
Rscript --vanilla -e "rmarkdown::render('03-add-cytoband-status-consensus.Rmd', clean = TRUE)"

# Run annotation step for consensus file
Rscript --vanilla 04-prepare-cn-file.R \
  --cnv_file ${scratch_dir}/consensus_seg_with_status.tsv \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "consensus_seg_annotated_cn" \
  --seg

# Define most focal units of recurrent CNVs
Rscript --vanilla -e "rmarkdown::render('05-define-most-focal-cn-units.Rmd', clean = TRUE)"

# Define the recurrent calls
Rscript --vanilla -e "rmarkdown::render('06-find-recurrent-calls.Rmd', clean = TRUE)"

libraryStrategies=("polya" "stranded")
chromosomesType=("autosomes" "x_and_y")
for strategy in ${libraryStrategies[@]}; do

  for chromosome_type in ${chromosomesType[@]}; do

    Rscript --vanilla rna-expression-validation.R \
      --annotated_cnv_file results/consensus_seg_annotated_cn_${chromosome_type}.tsv.gz \
      --expression_file ${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.${strategy}.rds \
      --independent_specimens_file $independent_specimens_file \
      --metadata $histologies_file \
      --goi_list $goi_file \
      --filename_lead "consensus_seg_annotated_cn"_${chromosome_type}_${strategy}
  done
done

# if we want to process the CNV data from the original callers
# (e.g., CNVkit, ControlFreeC)
if [ "$RUN_ORIGINAL" -gt "0" ]; then

  # Prep the CNVkit data
  Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"

  # Run annotation step for CNVkit
  Rscript --vanilla 04-prepare-cn-file.R \
    --cnv_file ${scratch_dir}/cnvkit_with_status.tsv \
    --gtf_file $gtf_file \
    --metadata $histologies_file \
    --filename_lead "cnvkit_annotated_cn" \
    --seg

  # Run annotation step for ControlFreeC
  Rscript --vanilla 04-prepare-cn-file.R \
    --cnv_file ${data_dir}/pbta-cnv-controlfreec.tsv.gz \
    --gtf_file $gtf_file \
    --metadata $histologies_file \
    --filename_lead "controlfreec_annotated_cn" \
    --controlfreec

  filenameLead=("cnvkit_annotated_cn" "controlfreec_annotated_cn")
  for filename in ${filenameLead[@]}; do
    for strategy in ${libraryStrategies[@]}; do
      for chromosome_type in ${chromosomesType[@]}; do
        Rscript --vanilla rna-expression-validation.R \
          --annotated_cnv_file results/${filename}_${chromosome_type}.tsv.gz \
          --expression_file ${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.${strategy}.rds \
          --independent_specimens_file $independent_specimens_file \
          --metadata $histologies_file \
          --goi_list $goi_file \
          --filename_lead ${filename}_${chromosome_type}_${strategy}
      done
    done
  done

fi
