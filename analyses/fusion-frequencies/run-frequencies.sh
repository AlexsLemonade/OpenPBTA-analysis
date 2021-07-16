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

# gather frequencies at FusionName and Fusion_Type level
Rscript 01-fusion-frequencies.R --fusion_file ../fusion_filtering/results/fusion-putative-oncogenic.tsv \
	--alt_id "FusionName,Fusion_Type" \
	--input_histologies ../../data/histologies.tsv \
	--primary_independence_list ../independent-samples/results/independent-specimens.rnaseq.primary.tsv \
	--relapse_independence_list ../independent-samples/results/independent-specimens.rnaseq.relapse.tsv \
	--output_filename "putative-oncogene-fusion-freq"

# gather frequencies at Fused Gene level
Rscript 01-fusion-frequencies.R --fusion_file ../fusion_filtering/results/fusion-putative-oncogenic.tsv \
        --alt_id "Gene_Symbol" \
        --input_histologies ../../data/histologies.tsv \
        --primary_independence_list ../independent-samples/results/independent-specimens.rnaseq.primary.tsv \
        --relapse_independence_list ../independent-samples/results/independent-specimens.rnaseq.relapse.tsv \
        --output_filename "putative-oncogene-fused-gene-freq"
        
# gather additional annotation for Fusion level         
Rscript 02-frequencies-annot.R --freq_file results/putative-oncogene-fusion-freq.tsv \
        --efo_mondo_file ../../data/efo-mondo-map.tsv  \
        --ensg_hugo_rmtl_file ../../data/ensg-hugo-rmtl-v1-mapping.tsv \
        --output_filename "putative-oncogene-fusion-freq-annotated"

# gather additional annotation for Fused Gene level           
Rscript 02-frequencies-annot.R --freq_file results/putative-oncogene-fused-gene-freq.tsv \
        --efo_mondo_file ../../data/efo-mondo-map.tsv  \
        --ensg_hugo_rmtl_file ../../data/ensg-hugo-rmtl-v1-mapping.tsv \
        --output_filename "putative-oncogene-fused-gene-freq-annotated"        

jq --compact-output '.[]' \
  results/putative-oncogene-fusion-freq-annotated.json \
  > results/putative-oncogene-fusion-freq-annotated.jsonl

jq --compact-output '.[]' \
  results/putative-oncogene-fused-gene-freq-annotated.json \
  > results/putative-oncogene-fused-gene-freq-annotated.jsonl

