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

# mkdir -p annotation_files
# wget --no-clobber --directory-prefix=annotation_files ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
# gunzip annotation_files/Homo_sapiens.GRCh38.84.gtf.gz

# per https://www.biostars.org/p/206342/
# awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' annotation_files/Homo_sapiens.GRCh38.84.gtf | gtf2bed - > annotation_files/Homo_sapiens.GRCh38.84.bed

Rscript --vanilla -e "rmarkdown::render('00-add-ploidy-cnvkit.Rmd', clean = TRUE)"
