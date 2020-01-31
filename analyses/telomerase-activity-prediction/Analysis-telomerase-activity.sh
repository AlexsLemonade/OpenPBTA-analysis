#enviroment settings
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mkdir plots

Rscript --vanilla 02-Comparing-Counts-versus-FPKM.R
Rscript --vanilla 03-Comparing-TERTexp-TERCexp-EXTENDScores.R
Rscript --vanilla 04-Comparing-Histology-versus-EXTENDScores.R