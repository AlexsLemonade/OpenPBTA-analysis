##########################################
#Purpose: Code to analyze and plot ssGSEA Results
#Author: Pichai Raman
#Date: 9/23/2019
##########################################


#Call libraries
library("tidyverse");

#Read in data
clinData <- read.delim("../../data/pbta-histologies.tsv");
geneSetExpMat <- readRDS("results/GeneSetExpressionMatrix.RDS")
