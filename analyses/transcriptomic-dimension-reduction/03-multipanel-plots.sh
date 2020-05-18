# Bethell and Taroni for CCDL 2019
# This generates multipanel dimension reduction plots for RSEM and kallisto
# data, with points colored either by RNA library or broad histology
#
# Usage: bash 03-multipanel-plots.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# loop over all RDS files that contain plot lists and save in the plots folder
FILES=plots/plot_data/*
for f in $FILES
do
  echo "Plotting $f ..."
  Rscript scripts/generate-multipanel-plot.R \
  --plot_rds $f \
  --plot_directory plots
done
