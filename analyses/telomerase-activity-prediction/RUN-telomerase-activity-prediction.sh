#enviroment settings
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

#generate collapsed data for count files 
Rscript ../analyses/collapse-rnaseq/01-summarize_matrices.R -i ../../data/pbta-gene-counts-rsem-expected_count.stranded.rds -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz -m ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds -t ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds
Rscript ../analyses/collapse-rnaseq/01-summarize_matrices.R -i ../../data/pbta-gene-counts-rsem-expected_count.polya.rds -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz -m ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds -t ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds


Rscript --vanilla 01-run-EXTEND.R --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_FPKM.txt
Rscript --vanilla 01-run-EXTEND.R --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output results/TelomeraseScores_PTBAPolya_FPKM.txt
Rscript --vanilla 01-run-EXTEND.R --input ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_counts.txt
Rscript --vanilla 01-run-EXTEND.R --input ../analyses/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds --output results/TelomeraseScores_PTBAPolya_counts.txt
#Rscript 01-run-EXTEND.R --input pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output TelomeraseScores_PTBAStranded_FPKM.txt