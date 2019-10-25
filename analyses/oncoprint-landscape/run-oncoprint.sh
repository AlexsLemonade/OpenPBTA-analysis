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

mkdir driver-lists
wget --directory-prefix=driver-lists -O driver-lists/brain-goi-list-long.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list-old.txt 
wget --directory-prefix=driver-lists -O driver-lists/brain-goi-list-short.txt https://github.com/marislab/create-pptc-pdx-oncoprints/raw/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data/brain-goi-list.txt

Rscript --vanilla 01-plot-oncoprint.R \
	--maf_file ../../data/pbta-snv-mutect2.vep.maf.gz \
	--fusion_file ../../data/pbta-fusion-arriba.tsv.gz 
