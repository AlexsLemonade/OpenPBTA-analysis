# Chante Bethell for CCDL 2019
# Run 01-prepare-cn-file.R
#
# Usage: bash run-prepare-cn.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript --vanilla 01-prepare-cn-file.R \
	--seg_file ../../data/pbta-cnv-cnvkit.seg.gz
