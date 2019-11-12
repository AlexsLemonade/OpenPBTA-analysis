# Chante Bethell for CCDL 2019
# Run 01-plot-oncoprint.R
#
# Usage: bash run-oncoprint.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mkdir -p driver-lists
wget --no-clobber --directory-prefix=driver-lists -O driver-lists/brain-goi-list-long.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list-old.txt
wget --no-clobber --directory-prefix=driver-lists -O driver-lists/brain-goi-list-short.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list.txt

# TODO: remove once the consensus files are included in the data download
# wget only the newer files
wget -N https://open-pbta.s3.amazonaws.com/data/snv-consensus/snv-consensus_11122019.zip
unzip -u snv-consensus_11122019.zip
# remove the mac specific directory included in the zip file
rm -r __MACOSX

Rscript --vanilla 01-plot-oncoprint.R \
	--maf_file ../../data/pbta-snv-mutect2.vep.maf.gz \
	--cnv_file ../../analyses/focal-cn-file-preparation/results/controlfreec_annotated_cn_autosomes.tsv.gz \
	--fusion_file ../../scratch/arriba.tsv \
	--goi_list driver-lists/brain-goi-list-long.txt \
	--png_name maf_oncoprint.png
