## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls ../../scratch/*.targetIntervals.bed`; do out=$(echo $i | sed 's/.bed/.Crossmapped_to_Gh38.bed/g'); CrossMap.py bed ../../scratch/hg19ToHg38.over.chain.gz $i $out; done

## This command sorts and merges all BED files in case there  are any overlaps in the BED files 
for i in `ls ../../scratch/*.targetIntervals.Crossmapped_to_Gh38.bed`; do out=$(echo $i | sed 's/.Crossmapped_to_Gh38.bed/.Gh38.bed/g');bedtools sort -i $i | bedtools merge  -i - > $out; done




