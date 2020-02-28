1) How did we get to this BED file from TCGA for each BAM file -
	a) Using this BAM file as an example (https://portal.gdc.cancer.gov/files/ffcec42c-999a-4b27-9e0e-f2f84f42acbe), use the case UUID to  query  for BED file(UUID also in the link)
	b) Use this curl statement to fetch target capture  information, we got the agilent BED file - 
		 curl 'https://api.gdc.cancer.gov/files/ffcec42c-999a-4b27-9e0e-f2f84f42acbe?fields=analysis.metadata.read_groups.target_capture_kit_name,analysis.metadata.read_groups.target_capture_kit_target_region&pretty=true'
	c) Results for curl statement from step above will get you to this BED file (https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed). Use wget  to download this BED file 
	d) Using liftover online, upload BED file and download the Gh38 BED(Removed additional contig chromosome lines using sed -i '/KI2/d' .bed) 


