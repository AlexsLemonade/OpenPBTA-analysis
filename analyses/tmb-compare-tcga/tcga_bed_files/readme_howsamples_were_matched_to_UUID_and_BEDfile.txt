1) MATCHING UP BED FILE NAMES AND UUID's
for i in `cat 160_caseids.txt`;  do echo "echo $i >> UUID_andBED.txt" >> curl.sh;echo "curl 'https://api.gdc.cancer.gov/files/"$i"?fields=analysis.metadata.read_groups.target_capture_kit_target_region&pretty=true' | grep "target_capture_kit_target_region" | sort | uniq | rev | cut -d "/" -f1 | rev >> UUID_andBED.txt" >> curl.sh; done

To run curl statements - 
sh curl.sh 

To move  BED file and UUIS to the same line
sed -i 's/\nwhole/\twhole/g'

2) OBTAINING BED FILES -
wget https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
wget https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/whole_exome_agilent_designed_120.targetIntervals.bed
wget https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed


3) Combining BED_UUID with UUID_sample with this file - 
python3 match_files.py

4) Used liftover online to convert all three BED files 
sed -i.bak 's/^/chr/' whole_exome_agilent_designed_120.targetIntervals.bed 
sed -i.bak 's/^/chr/'  whole_exome_agilent_plus_tcga_6k.targetIntervals.bed
bedtools sort -i whole_exome_agilent_designed_120.targetIntervals.bed > whole_exome_agilent_designed_120.targetIntervals.sorted.bed
bedtools sort -i whole_exome_agilent_plus_tcga_6k.targetIntervals.bed > whole_exome_agilent_plus_tcga_6k.targetIntervals.sorted.bed
sed -i.bak '/KI2/d' whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed
sed -i.bak '/KI2/d' whole_exome_agilent_designed_120.targetIntervals.Gh38.bed





