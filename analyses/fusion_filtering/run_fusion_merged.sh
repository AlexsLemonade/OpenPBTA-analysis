#!/bin/bash

# K S Gaonkar

# Run fusion_filtering

set -e
set -o pipefail

# fusion files before and after standardization
arriba_file="data/pbta-fusion-arriba.tsv.gz"
starfusion_file="data/pbta-fusion-starfusion.tsv.gz"
standard_arriba_file="scratch/arriba.tsv"
standard_starfusion_file="scratch/starfusion.tsv"

# general filtering parameters
artifact_filter="GTEx_recurrent_STARF2019|HGNC_GENEFAM|DGD_PARALOGS|Normal|BodyMap|ConjoinG"
reading_frame_filter="in-frame|frameshift|other"

# relevant gene expression files
polya_expression_file="data/pbta-gene-expression-rsem-fpkm.polya.rds"
stranded_expression_file="data/pbta-gene-expression-rsem-fpkm.stranded.rds"

# directory that holds all of the reference files
references_directory="analyses/fusion_filtering/references/"
normal_expression_file="$references_directory/Brain_FPKM_hg38_matrix.txt.zip"

# Run Fusion standardization
Rscript analyses/fusion_filtering/01-fusion-standardization.R --fusionfile $arriba_file \
                                                              --caller "arriba" \
                                                              --outputfile $standard_arriba_file
Rscript analyses/fusion_filtering/01-fusion-standardization.R --fusionfile $starfusion_file \
                                                              --caller "starfusion" \
                                                              --outputfile $standard_starfusion_file

# Run Fusion general filtering
Rscript analyses/fusion_filtering/02-fusion-filtering.R --standardFusionFileSTARFusion $standard_starfusion_file  \
							--standardFusionFileArriba $standard_arriba_file \
                                                        --expressionMatrix $polya_expression_file \
                                                        --artifactFilter $artifact_filter  \
                                                        --readingFrameFilter $reading_frame_filter \
                                                        --referenceFolder $references_directory \
                                                        --outputfile scratch/standardFusionPolyaExp \
                                                        --readthroughFilter

Rscript analyses/fusion_filtering/02-fusion-filtering.R --standardFusionFileArriba $standard_arriba_file \
							--standardFusionFileSTARFusion $standard_starfusion_file  \
                                                        --expressionMatrix $stranded_expression_file \
                                                        --artifactFilter $artifact_filter \
                                                        --readingFrameFilter $reading_frame_filter \
                                                        --referenceFolder $references_directory \
                                                        --outputfile scratch/standardFusionStrandedExp \
                                                        --readthroughFilter


# Fusion zscore annotation for filtered fusion
Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R --standardFusionCalls scratch/standardFusionPolyaExp_QC_expression_filtered_annotated.RDS \
                                                            --expressionMatrix $polya_expression_file \
                                                            --normalExpressionMatrix $normal_expression_file \
                                                            --outputfile scratch/standardFusionPolyaExp_QC_expression

Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R --standardFusionCalls scratch/standardFusionStrandedExp_QC_expression_filtered_annotated.RDS \
                                                            --expressionMatrix $stranded_expression_file \
                                                            --normalExpressionMatrix $normal_expression_file \
                                                            --outputfile scratch/standardFusionStrandedExp_QC_expression

