#!/bin/bash

# K S Gaonkar

# Run fusion_filtering

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up paths to data files consumed by analysis, and path to result output
data_path="../../data"
scratch_path="../../scratch"
references_path="references"
results_path="results/"


# fusion files before and after standardization
arriba_file="${data_path}/pbta-fusion-arriba.tsv.gz"
starfusion_file="${data_path}/pbta-fusion-starfusion.tsv.gz"
standard_arriba_file="${scratch_path}/arriba.tsv"
standard_starfusion_file="${scratch_path}/starfusion.tsv"

# general filtering parameters
artifact_filter="GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG"
reading_frame_filter="in-frame|frameshift|other"

# relevant gene expression files
polya_expression_file="${data_path}/pbta-gene-expression-rsem-fpkm.polya.rds"
stranded_expression_file="${data_path}/pbta-gene-expression-rsem-fpkm.stranded.rds"

# reference expression file
normal_expression_file="${references_path}/Brain_FPKM_hg38_matrix.txt.zip"

# metadata files
histologies_file="${data_path}/pbta-histologies.tsv"
independent_samples_file="${data_path}/independent-specimens.wgswxs.primary-plus.tsv"

# data release files to use for recurrent fusion/fused genes detection

putative_oncogenic_fusion="${data_path}/pbta-fusion-putative-oncogenic.tsv"


# Run Fusion standardization for arriba caller
Rscript 01-fusion-standardization.R --fusionfile $arriba_file \
                                    --caller "arriba" \
                                    --outputfile $standard_arriba_file
                                    
                                    
# Run Fusion standardization for starfusion caller
Rscript 01-fusion-standardization.R --fusionfile $starfusion_file \
                                    --caller "starfusion" \
                                    --outputfile $standard_starfusion_file

# Run Fusion general filtering for polya
Rscript 02-fusion-filtering.R --standardFusionFiles $standard_starfusion_file,$standard_arriba_file  \
                              --expressionMatrix $polya_expression_file \
                              --artifactFilter $artifact_filter  \
                              --readingFrameFilter $reading_frame_filter \
                              --referenceFolder $references_path \
                              --outputfile "${scratch_path}/standardFusionPolyaExp" \
                              --readthroughFilter
                              
                              
# Run Fusion general filtering for stranded
Rscript 02-fusion-filtering.R --standardFusionFiles $standard_arriba_file,$standard_starfusion_file \
                              --expressionMatrix $stranded_expression_file \
                              --artifactFilter $artifact_filter \
                              --readingFrameFilter $reading_frame_filter \
                              --referenceFolder $references_path \
                              --outputfile "${scratch_path}/standardFusionStrandedExp" \
                              --readthroughFilter


# Fusion zscore annotation for filtered fusion for polya
Rscript 03-Calc-zscore-annotate.R --standardFusionCalls "${scratch_path}/standardFusionPolyaExp_QC_expression_filtered_annotated.RDS" \
                                  --expressionMatrix $polya_expression_file \
                                  --normalExpressionMatrix $normal_expression_file \
                                  --outputfile "${scratch_path}/standardFusionPolyaExp_QC_expression"

# Fusion zscore annotation for filtered fusion for stranded
Rscript 03-Calc-zscore-annotate.R --standardFusionCalls "${scratch_path}/standardFusionStrandedExp_QC_expression_filtered_annotated.RDS" \
                                  --expressionMatrix $stranded_expression_file \
                                  --normalExpressionMatrix $normal_expression_file \
                                  --outputfile "${scratch_path}/standardFusionStrandedExp_QC_expression"

# Project specific filtering
Rscript -e "rmarkdown::render('04-project-specific-filtering.Rmd')"

# QC filter putative oncogene found in more than 4 histologies
Rscript -e "rmarkdown::render('05-QC_putative_onco_fusion_dustribution.Rmd')"

# Recurrent fusion/fused genes
Rscript 06-recurrent-fusions-per-histology.R --standardFusionCalls $putative_oncogenic_fusion \
                                             --clinicalFile $histologies_file \
                                             --outputfolder $results_path \
                                             --independentSpecimensFile $independent_samples_file
