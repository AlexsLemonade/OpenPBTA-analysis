#!/bin/bash

#download fusion filtering git repo to access the scripts and refernce files
cd ../../scratch
git clone https://github.com/d3b-center/fusion_filtering_pipeline.git


#get back to analysis folder
cd ../../analyses/fusion_filtering/

#copy Rscripts from git repo
cp ../../scratch/fusion_filtering_pipeline/R/* ./

#driver list:filter for biologically meaningful gene list
#filtered list: fusions called by both callers, unique to broad histology, in atleast 2 samples per broad histology
Rscript 02-fusion-filtering.R

#filt for low expressed gene1 and gene2 in driver list
Rscript 03-filter-low-exp-driver.R

#filter out low expressed gene1 and gene2 in filtered list
Rscript 04-filter-low-exp-prioritised.R

#annotated filtered fusion with TSG/Oncogenic/Kinase/TF and TCGA calls
Rscript 05-annotate-filtered-fusion.R

# Plotting recurrent driver list fusion and fused genes
Rscript 06-plot-recurrent-fusions.R

