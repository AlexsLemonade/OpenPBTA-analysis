# This script downloads the chain file from UCSC, converts the  BED files from hg19 to Gh38 ans sorts and mereges (in case there are  any overlaps in the BED regions)


## Downloading chain.gz file for CrossMap command 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz  -O ../../scratch/hg19ToHg38.over.chain.gz

## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls ../../scratch/*.targetIntervals.bed`; do 
  out=$(echo $i | sed 's/.bed/.Crossmapped_to_Gh38.bed/g')
  CrossMap.py bed ../../scratch/hg19ToHg38.over.chain.gz $i $out
done


## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls ../../scratch/*.targetIntervals.Gh38.bed`; do 
  bedtools sort -i $i \
  | bedtools merge \
  > results/$(basename $i)
done



