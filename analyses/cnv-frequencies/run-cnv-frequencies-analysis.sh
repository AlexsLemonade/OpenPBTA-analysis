#!/bin/bash
# Eric Wafula
# Pediatric OpenTargets 2021

# Purpose: Run a CNV consensus gene-level frequecies anaysis for PediatricOpenTargets

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# All file paths in this bash script are based on root directory of this Git repository,
# therefore the script should always be run from the root directory of OPenPedCan-analysis
# Adapted from OpenPBTA analysis modules


# Histology file path 
histology_file=data/histologies.tsv

# CNV consensus file path
cnv_file=data/consensus_seg_annotated_cn_autosomes.tsv

# Independent primary tumor samples file path
primary_tumors=analyses/independent-samples/results/independent-specimens.wgs.primary.tsv

# Independent relapse tumor samples file path
relapse_tumors=analyses/independent-samples/results/independent-specimens.wgs.relapse.tsv

# Disease to EFO/MONDO mapping file path
efo_mondo=data/efo-mondo-map.tsv

# Ensembl to RMTL mapping file path
ensg_rmtl=data/ensg-hugo-rmtl-v1-mapping.tsv


python3 analyses/cnv-frequencies/01-cnv-frequencies.py \
 $histology_file $cnv_file $primary_tumors $relapse_tumors $efo_mondo $ensg_rmtl
