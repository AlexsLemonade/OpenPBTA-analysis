#!/bin/bash

# K S Gaonkar

# Run fusion_filtering
# Takes one environment variable, `OPENPBTA_BASE_SUBTYPING`, if value is 1 then
# uses pbta-histologies-base.tsv for subtyping if value is 0 runs all modules (Default)
# with pbta-histologies.tsv

set -e
set -o pipefail

RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

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
arriba_file="${data_path}/fusion-arriba.tsv.gz"
starfusion_file="${data_path}/fusion-starfusion.tsv.gz"
standard_arriba_file="${scratch_path}/arriba.tsv"
standard_starfusion_file="${scratch_path}/starfusion.tsv"

# general filtering parameters
artifact_filter="GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap"
reading_frame_filter="in-frame|frameshift|other"
spanningFragCountFilter=100


# relevant gene expression files
rna_expression_file="${data_path}/gene-expression-rsem-tpm-collapsed.rds"

# reference expression file
normal_expression_adrenal_gland="${references_path}/gtex_adrenal_gland_TPM_hg38.rds"
normal_expression_brain="${references_path}/gtex_brain_TPM_hg38.rds"

# metadata files
if [[ RUN_FOR_SUBTYPING -eq "0" ]]
then
   histologies_file="${data_path}/histologies.tsv" 
   independent_samples_file="${data_path}/independent-specimens.wgswxspanel.primary-plus.tsv"
else 
   histologies_file="${data_path}/histologies-base.tsv"  
   independent_samples_file="../independent-samples/results/independent-specimens.wgswxspanel.primary-plus.tsv" 
fi


# data release files to use for recurrent fusion/fused genes detection

putative_oncogenic_fusion="${results_path}/pbta-fusion-putative-oncogenic.tsv"

# Run filtering code to get the reference file
Rscript 00-normal-matrix-generation.R  --expressionMatrix $rna_expression_file \
                                       --clinicalFile $histologies_file \
                                       --specimenType "Adrenal Gland" \
                                       --outputFile $normal_expression_adrenal_gland


Rscript 00-normal-matrix-generation.R  --expressionMatrix $rna_expression_file \
                                       --clinicalFile $histologies_file \
                                       --specimenType "Brain" \
                                       --outputFile $normal_expression_brain

# Run Fusion standardization for arriba caller
Rscript 01-fusion-standardization.R --fusionfile $arriba_file \
                                    --caller "arriba" \
                                    --outputFile $standard_arriba_file
                                    
                                    
# Run Fusion standardization for starfusion caller
Rscript 01-fusion-standardization.R --fusionfile $starfusion_file \
                                    --caller "starfusion" \
                                    --outputFile $standard_starfusion_file

# Run Fusion general filtering for combined expression file
Rscript 02-fusion-filtering.R --standardFusionFiles $standard_starfusion_file,$standard_arriba_file  \
                              --expressionMatrix $rna_expression_file \
                              --clinicalFile $histologies_file \
                              --cohortInterest "PBTA,GMKF" \
                              --artifactFilter $artifact_filter  \
                              --spanningFragCountFilter $spanningFragCountFilter \
                              --readingFrameFilter $reading_frame_filter \
                              --referenceFolder $references_path \
                              --outputFile "${scratch_path}/standardFusionExp" \
                              --readthroughFilter


# Fusion zscore annotation for filtered fusion for the combined RNA expression file
Rscript 03-Calc-zscore-annotate.R --standardFusionCalls "${scratch_path}/standardFusionExp_QC_expression_filtered_annotated.RDS" \
                                  --expressionMatrix $rna_expression_file \
                                  --clinicalFile $histologies_file \
                                  --cohortInterest "PBTA,GMKF" \
                                  --normalExpressionMatrix $normal_expression_brain,$normal_expression_adrenal_gland \
                                  --outputFile "${scratch_path}/standardFusionExp_QC_expression"

# Project specific filtering
Rscript -e "rmarkdown::render('04-project-specific-filtering.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# QC filter putative oncogene found in more than 4 histologies
Rscript -e "rmarkdown::render('05-QC_putative_onco_fusion_distribution.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# Recurrent fusion/fused genes
Rscript 06-recurrent-fusions-per-cancer-group.R --standardFusionCalls $putative_oncogenic_fusion \
                                                --clinicalFile $histologies_file \
                                                --cohortInterest "PBTA,GMKF" \
                                                --outputfolder $results_path \
                                                --independentSpecimensFile $independent_samples_file

