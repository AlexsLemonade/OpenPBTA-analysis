#!/bin/bash


# This script captures the steps to generate the bad_chromosomal_seg_updated_merged.txt 

# There are two components to this file the "IGLL regions, centromeric and telomeric regions" and the "segmental duplication regions"
# 1) The IGLL regions, centromeric and telomeric regions are generated from the practice described by Kai Wang at his PennCNV website http://penncnv.openbioinformatics.org/en/latest/misc/faq/
# 2) The segmental duplication are downloaded from UCSC genome browser. The segmental duplication with 95% identity was downloaded and merged. 
# 3) In the end, everything was put together into one file, sorted and merged. 



cat 
