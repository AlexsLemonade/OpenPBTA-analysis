samtools depth -b whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.liftover_Gh38.bed --reference /sbgenomics/project-files/Homo_sapiens_assembly38.fasta /sbgenomics/project-files/C494.TCGA-HT-7684-01A-11D-2253-08.1_gdc_realn.bam  > depth_with_whole_exome_agilent_1.1_refseq_plus_3_boosters_TCGA-HT-7684-01A-11D-2253-08

samtools depth -b whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.liftover_Gh38.bed --reference /sbgenomics/project-files/Homo_sapiens_assembly38.fasta /sbgenomics/project-files/C1663.TCGA-SP-A6QF-01A-12D-A35I-08.1_gdc_realn.bam > depth_with_whole_exome_agilent_1.1_refseq_plus_3_boosters_TCGA-SP-A6QF-01A-12D-A35I-0

samtools depth -b /sbgenomics/project-files/gencode.v19.basic.exome_withchr_Gh38_sorted_merged.bed --reference /sbgenomics/project-files/Homo_sapiens_assembly38.fasta /sbgenomics/project-files/C494.TCGA-HT-7684-01A-11D-2253-08.1_gdc_realn.bam > depth_with_gencode_Gh38_sorted_merged_TCGA-HT-7684-01A-11D-2253-08

samtools depth -b /sbgenomics/project-files/gencode.v19.basic.exome_withchr_Gh38_sorted_merged.bed --reference /sbgenomics/project-files/Homo_sapiens_assembly38.fasta /sbgenomics/project-files/C1663.TCGA-SP-A6QF-01A-12D-A35I-08.1_gdc_realn.bam > depth_with_gencode_Gh38_sorted_merged_TCGA-SP-A6QF-01A-12D-A35I-08






Numbers in BED file - 
awk '{if($3>20) print}' depth_with_whole_exome_agilent_1.1_refseq_plus_3_boosters_TCGA-SP-A6QF-01A-12D-A35I-0 | wc -l
27350549
awk '{if($3>20) print}' depth_with_whole_exome_agilent_1.1_refseq_plus_3_boosters_TCGA-HT-7684-01A-11D-2253-08 | wc -l
27755836

Count BED regions - 
awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}'  /sbgenomics/project-files/gencode.v19.basic.exome_withchr_Gh38_sorted_merged.bed
94468800
awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.liftover_Gh38.bed
32976776
