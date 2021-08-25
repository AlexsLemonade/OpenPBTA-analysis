#!/bin/bash
# PediatricOpenTargets 2021
# Yuanchao Zhang
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# copied from https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/scripts/run_in_ci.sh
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up paths to data files consumed by analysis
data_path="../../data"

# Independent sample lists needed 
isl_primary_each="${data_path}/independent-specimens.rnaseq.primary.eachcohort.tsv"
isl_relapse_each="${data_path}/independent-specimens.rnaseq.relapse.eachcohort.tsv"
isl_primary_all="${data_path}/independent-specimens.rnaseq.primary.tsv"
isl_relapse_all="${data_path}/independent-specimens.rnaseq.relapse.tsv"

# Filtered Fusion file 
fusion_file="${data_path}/fusion-putative-oncogenic.tsv"

# Histology file
histology_file="${data_path}/histologies.tsv"

# gather frequencies at FusionName and Fusion_Type level
Rscript 01-fusion-frequencies.R --fusion_file $fusion_file \
	--alt_id "FusionName,Fusion_Type" \
	--input_histologies $histology_file \
	--primary_independence_all $isl_primary_all \
	--relapse_independence_all $isl_relapse_all \
	--primary_independence_each $isl_primary_each \
	--relapse_independence_each $isl_relapse_each \
	--output_filename "putative-oncogene-fusion-freq" 

jq --compact-output '.[]' \
  results/putative-oncogene-fusion-freq.json \
  > results/putative-oncogene-fusion-freq.jsonl

rm results/putative-oncogene-fusion-freq.json

# gather frequencies at Fused Gene level
Rscript 01-fusion-frequencies.R --fusion_file $fusion_file \
        --alt_id "Gene_Symbol" \
        --input_histologies $histology_file \
        --primary_independence_all $isl_primary_all \
        --relapse_independence_all $isl_relapse_all \
      	--primary_independence_each $isl_primary_each \
      	--relapse_independence_each $isl_relapse_each \
        --output_filename "putative-oncogene-fused-gene-freq" 

jq --compact-output '.[]' \
  results/putative-oncogene-fused-gene-freq.json \
  > results/putative-oncogene-fused-gene-freq.jsonl

rm results/putative-oncogene-fused-gene-freq.json
gzip results/putative-oncogene*        
