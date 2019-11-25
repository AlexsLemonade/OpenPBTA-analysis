#!/bin/bash
# This script runs the sex-prediction-from-RNASeq analysis
# Author's Name Bill Amadio 2019

set -e
set -o pipefail

Rscript --vanilla 01-clean_split_data.R --expression ../../data/pbta-gene-expression-kallisto.stranded.rds --metadata ../../data/pbta-histologies.tsv --output_directory processed_data --filename_lead kallisto_stranded

Rscript --vanilla 02-train_elasticnet.R --input_directory processed_data --output_directory models --filename_lead kallisto_stranded

Rscript --vanilla 03-evaluate_model.R --test_set_input_directory processed_data --model_input_directory models --output_directory results --filename_lead kallisto_stranded

