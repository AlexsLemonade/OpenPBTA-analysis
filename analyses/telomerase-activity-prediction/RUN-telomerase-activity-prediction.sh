#enviroment settings
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# If this is set to 1, only run code necessary for generating manuscript figures: 01-run-EXTEND.R on stranded, as well as the thresholded analysis.
RUN_FOR_FIGURES=${OPENPBTA_FOR_FIGURES:-0}


mkdir -p results
mkdir -p plots

if [ "${RUN_FOR_FIGURES}" == "1" ]; then

  #generate telomerase activities using gene expression data from collapse RNA seq data files
  Rscript --vanilla 01-run-EXTEND.R --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_FPKM.txt

else

  #generate collapsed data for count files
  Rscript ../collapse-rnaseq/01-summarize_matrices.R -i ../../data/pbta-gene-counts-rsem-expected_count.stranded.rds -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz -m ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds -t ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds
  Rscript ../collapse-rnaseq/01-summarize_matrices.R -i ../../data/pbta-gene-counts-rsem-expected_count.polya.rds -g ../../data/gencode.v27.primary_assembly.annotation.gtf.gz -m ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds -t ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds

  #generate telomerase activities using gene expression data from collapse RNA seq data files
  Rscript --vanilla 01-run-EXTEND.R --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_FPKM.txt
  Rscript --vanilla 01-run-EXTEND.R --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output results/TelomeraseScores_PTBAPolya_FPKM.txt
  Rscript --vanilla 01-run-EXTEND.R --input ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_counts.txt
  Rscript --vanilla 01-run-EXTEND.R --input ../collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds --output results/TelomeraseScores_PTBAPolya_counts.txt

  #Compare Telomerase scores for PolyA and Stranded data
  echo "Comparing telomerase scores for polyA and stranded data..."
  Rscript --vanilla 02-Comparing-Counts-versus-FPKM.R --output plots/PTBA_GE_Score_AllScatter.pdf

  #Compare Telomerase scores with TERT and TERC gene expressions
  echo "Comparing EXTEND scores with TERT and TERC expression..."
  Rscript --vanilla 03-Comparing-TERTexp-TERCexp-EXTENDScores.R --output plots/PTBA_GE_TM_ScatterComp.pdf

  #Distribution of telomerase scores across various histologies of brain tumors.
  echo "Plotting distribution of EXTEND scores in histologies..."
  Rscript --vanilla 04-Comparing-Histology-versus-EXTENDScores.R --output plots/PBTA_StrandedHistology.pdf

  #Distribution of telomerase scores across different molecular subtypes of medulloblastoma tumors.
  echo "Plotting distribution of EXTEND scores in MB subtypes..."
  Rscript --vanilla 05-Comparing-MolecularSubtypes-EXTENDScores.R

  # Distribution of telomerase scores in TERTp-positive vs. TERTp-negative samples
  echo "Plotting distribution of EXTEND scores in samples with or without TERTp..."
  Rscript --vanilla 06-TERTp-EXTENDScores.R
fi

# Regardless of which run option, we would like to analyze results at a tumor purity threshold:
echo "Exploring results at a tumor purity threshold..."

# Generate EXTEND scores
Rscript --vanilla 01-run-EXTEND.R --apply_tumor_purity_threshold --input ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output results/TelomeraseScores_PTBAStranded_FPKM_thresholded.txt

# Explore thresholded results
Rscript --vanilla -e "rmarkdown::render('07_EXTEND-at-threshold.Rmd')"
