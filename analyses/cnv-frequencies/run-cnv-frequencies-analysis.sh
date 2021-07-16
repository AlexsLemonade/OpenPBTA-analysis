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
autosomes_cnv_file=data/consensus_seg_annotated_cn_autosomes.tsv.gz
allosomes_cnv_file=data/consensus_seg_annotated_cn_x_and_y.tsv.gz

# Independent primary tumor samples file path
primary_tumors=analyses/independent-samples/results/independent-specimens.wgs.primary.tsv

# Independent relapse tumor samples file path
relapse_tumors=analyses/independent-samples/results/independent-specimens.wgs.relapse.tsv

# Disease to EFO/MONDO mapping file path
efo_mondo=data/efo-mondo-map.tsv

# Ensembl to RMTL mapping file path
ensg_rmtl=data/ensg-hugo-rmtl-v1-mapping.tsv

# Hugo gene symbols to OncoKB categories mapping file path
oncokb=analyses/cnv-frequencies/input/OncoKB_Oncogene-TSG_genes.tsv

####### compute autosomes CNV frequencies #############
python3 analyses/cnv-frequencies/01-cnv-frequencies.py \
	$histology_file $autosomes_cnv_file $primary_tumors $relapse_tumors $oncokb $efo_mondo $ensg_rmtl

####### compute aullosomes CNV frequencies #############
python3 analyses/cnv-frequencies/01-cnv-frequencies.py \
	$histology_file $allosomes_cnv_file $primary_tumors $relapse_tumors $oncokb $efo_mondo $ensg_rmtl

####### compress the output files ######################
gzip analyses/cnv-frequencies/results/consensus_seg_annotated_cn_*
