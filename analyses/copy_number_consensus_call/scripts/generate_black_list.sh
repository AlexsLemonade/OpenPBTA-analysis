#!/bin/bash


# This script captures the steps to generate the "blacklist" regions where CNV calls are not likely to be reliable,
# and will be excluded from downstream analysis.
#
# 1) The Immunogobulin (IG) regions, centromeric and telomeric regions are generated from
# the practice described by Kai Wang at his PennCNV website http://penncnv.openbioinformatics.org/en/latest/misc/faq/
#  a) The IG regions are defined by the following coordinates described by Michael Xie (using hg18 coordinates):
#          chr2:88935000-89418000 IgKappa
#          chr6:29775000-33225000 HLA*
#          chr7:141636000-142225000 TCRbeta
#          chr14:21214600-22095500 TCRalpha
#          chr14:105046000-106368585 IgHeavy
#          chr22:20675000-21620000 IgLambda
#      These were converted to hg38 coordinates to generate the file `ref/immunoglobulin-regions.bed`
#      These regions are slightly different from those defined on Kai Wang's website, but are mostly in agreement.
#  b) Centromeric and heterchromtic regions are extracted from the following file at UCSC
#     http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz

curl -O  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
gunzip -c cytoBand.txt.gz | awk -v OFS='\t' '{if($5=="acen"){print}}' > ref/centromeres.bed
gunzip -c cytoBand.txt.gz | awk -v OFS='\t' '{if($5=="gvar"){print}}' > ref/heterochromatin.bed


# 2) The segmental duplication are downloaded from UCSC genome browser. Segmental duplication with 95% or greater identity are extracted and merged.
#    The track used for this is here https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=783211655_qnHWS1w9VWkUebtWb7482jKHTMF8&c=chr1&g=genomicSuperDups
#    http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

# segdups
curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
gunzip -c genomicSuperDups.txt.gz | awk -v OFS='\t' '{if($27>0.95){print $2, $3, $4, $27}}' |bedtools merge > genomicSuperDups95.bed


# After having the IGLL_telo_centromeric_region.txt from the method described by Kai Wang, the following were performed to form the black list
# Sort the segdup file and merge it.
bedtools merge > segdup_95_merged


# Put the IGLL_telo_centromeric_region.txt and the segdup file together,

# take only the first 3 columns,
# sort them and merge any overlaping segment
cat IGLL_telo_centromeric_region.txt segdup_95_merged | awk -v OFS='\t' '{print$1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge > bad_chromosomal_seg_updated_merged.bed
