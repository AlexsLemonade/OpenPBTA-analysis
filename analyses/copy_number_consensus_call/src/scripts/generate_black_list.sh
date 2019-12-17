#!/bin/bash


# This script captures the steps to generate the bad_chromosomal_seg_updated_merged.txt 

# There are two components to this file the "IGLL regions, centromeric and telomeric regions" and the "segmental duplication regions"
# 1) The IGLL regions, centromeric and telomeric regions are generated from the practice described by Kai Wang at his PennCNV website http://penncnv.openbioinformatics.org/en/latest/misc/faq/
# 2) The segmental duplication are downloaded from UCSC genome browser. The segmental duplication with 95% identity was downloaded and merged. 
# 3) In the end, everything was put together into one file, sorted and merged. 



# After having the IGLL_telo_centromeric_region.txt from the method described by Kai Wang, the following were performed to form the black list
# Sort the segdup file and merge it. 
sort -k1,1 -k2,2n segdup_95 | bedtools merge > segdup_95_merged


# Put the IGLL_telo_centromeric_region.txt and the segdup file together,
# take only the first 3 columns,
# sort them and merge any overlaping segment
cat IGLL_telo_centromeric_region.txt segdup_95_merged | awk -v OFS='\t' '{print$1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge > bad_chromosomal_seg_updated_merged.txt
