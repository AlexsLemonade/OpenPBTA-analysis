# Chante Bethell for CCDL 2019
# Run 01-plot-oncoprint.R
#
# Usage: bash run-oncoprint.sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

#### Download the driver lists from the PPTC repository

mkdir -p driver-lists
wget -N --directory-prefix=driver-lists -O driver-lists/brain-goi-list-long.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list-old.txt
wget -N --directory-prefix=driver-lists -O driver-lists/brain-goi-list-short.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list.txt

#### Download consensus mutation files

maf_consensus=../../data/pbta-snv-consensus-mutation.maf.tsv.gz
controlfreec_autosomes=../focal-cn-file-preparation/results/controlfreec_annotated_cn_autosomes.tsv.bz2
fusion_file=../../data/pbta-fusion-putative-oncogenic.tsv
histologies_file=../../data/pbta-histologies.tsv
intermediate_directory=../../scratch/oncoprint_files
primary_filename="all_participants_primary_only"
primaryplus_filename="all_participants_primary-plus"
genes_list=driver-lists/brain-goi-list-long.txt

#### Primary only oncoprint

Rscript --vanilla 00-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --cnv_file ${controlfreec_autosomes} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primary_filename} \
  --independent_specimens ../../data/independent-specimens.wgswxs.primary.tsv

Rscript --vanilla 01-plot-oncoprint.R \
  --maf_file ${intermediate_directory}/${primary_filename}_maf.tsv \
  --cnv_file ${intermediate_directory}/${primary_filename}_cnv.tsv \
  --fusion_file ${intermediate_directory}/${primary_filename}_fusions.tsv \
  --metadata_file ${histologies_file} \
  --png_name ${primary_filename}_oncoprint.png

# Genes of interest only version of oncoprint
Rscript --vanilla 01-plot-oncoprint.R \
  --maf_file ${intermediate_directory}/${primary_filename}_maf.tsv \
  --cnv_file ${intermediate_directory}/${primary_filename}_cnv.tsv \
  --fusion_file ${intermediate_directory}/${primary_filename}_fusions.tsv \
  --metadata_file ${histologies_file} \
  --goi_list ${genes_list} \
  --png_name ${primary_filename}_goi_oncoprint.png

#### Primary plus samples oncoprint

Rscript --vanilla 00-map-to-sample_id.R \
  --maf_file ${maf_consensus} \
  --cnv_file ${controlfreec_autosomes} \
  --fusion_file ${fusion_file} \
  --metadata_file ${histologies_file} \
  --output_directory ${intermediate_directory} \
  --filename_lead ${primaryplus_filename} \
  --independent_specimens ../../data/independent-specimens.wgswxs.primary-plus.tsv

Rscript --vanilla 01-plot-oncoprint.R \
  --maf_file ${intermediate_directory}/${primaryplus_filename}_maf.tsv \
  --cnv_file ${intermediate_directory}/${primaryplus_filename}_cnv.tsv \
  --fusion_file ${intermediate_directory}/${primaryplus_filename}_fusions.tsv \
  --metadata_file ${histologies_file} \
  --png_name ${primaryplus_filename}_oncoprint.png

# Genes of interest only version of oncoprint
Rscript --vanilla 01-plot-oncoprint.R \
  --maf_file ${intermediate_directory}/${primaryplus_filename}_maf.tsv \
  --cnv_file ${intermediate_directory}/${primaryplus_filename}_cnv.tsv \
  --fusion_file ${intermediate_directory}/${primaryplus_filename}_fusions.tsv \
  --metadata_file ${histologies_file} \
  --goi_list ${genes_list} \
  --png_name ${primaryplus_filename}_goi_oncoprint.png
