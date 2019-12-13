## Data file descriptions

This document contains information about all data files associated with this project. Each file should have the following association information:

+ **File type** should be one of..
	+ *Reference file*: Obtained from an external source/database. When known, the obtained data and a link to the external source should be included.
	+ *Modified reference file*: Obtained from an external source/database but modified for OpenPBTA use. 	
	+ *External data file*: Data directly obtained from the cancer samples databases that be. When known, the specific database and download date should be included.
	+ *Analysis file*: Any file created by a script in `analyses/*`. 
+ **Associated analyses**
	+ A relative link to the specific analysis from which the file was generated.
		+ For any files which are generally applicable to many/most/all analyses, please write *Universal* in this field
	+ When applicable, a link to the specific *script* that produced (or modified, for *Modified reference file* types) the data
+ **File description**
	+ A *brief* one sentence description of what the file contains (e.g., bed files contain coordinates for features XYZ).



### current release (release-v11-20191126)

| **File name** |  **File Type** | **Associated analysis** | **File Description** |
|---------------|----------------|------------------------|-----------------------|
|`GRCh38.primary_assembly.genome.fa.gz` | | |  |
|`StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed` | | |
|`StrexomeLite_hg38_liftover_100bp_padded.bed`| | |
|`WGS.hg38.lancet.300bp_padded.bed` | | |
|`WGS.hg38.lancet.unpadded.bed` | | |
|`WGS.hg38.mutect2.unpadded.bed` | | |
|`WGS.hg38.strelka2.unpadded.bed` | ||
|`WGS.hg38.vardict.100bp_padded.bed` | ||
|`WXS.hg38.100bp_padded.bed` | ||
|`gencode.v27.primary_assembly.annotation.gtf.gz` | ||
|`independent-specimens.wgs.primary-plus.tsv` | ||
|`independent-specimens.wgs.primary.tsv` | ||
|`independent-specimens.wgswxs.primary-plus.tsv` | ||
|`independent-specimens.wgswxs.primary.tsv` | ||
|`pbta-cnv-cnvkit.seg.gz` | ||
|`pbta-cnv-controlfreec.tsv.gz` | ||
|`pbta-fusion-arriba.tsv.gz` | ||
|`pbta-fusion-putative-oncogenic.tsv` | ||
|`pbta-fusion-starfusion.tsv.gz` | ||
|`pbta-gene-counts-rsem-expected_count.polya.rds` | ||
|`pbta-gene-counts-rsem-expected_count.stranded.rds` || |
|`pbta-gene-expression-kallisto.polya.rds` | ||
|`pbta-gene-expression-kallisto.stranded.rds` || |
|`pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` | ||
|`pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` || |
|`pbta-gene-expression-rsem-fpkm.polya.rds` | ||
|`pbta-gene-expression-rsem-fpkm.stranded.rds` || |
|`pbta-gene-expression-rsem-tpm.polya.rds` | ||
|`pbta-gene-expression-rsem-tpm.stranded.rds` || |
|`pbta-histologies.tsv` | ||
|`pbta-isoform-counts-rsem-expected_count.polya.rds` | ||
|`pbta-isoform-counts-rsem-expected_count.stranded.rds` || |
|`pbta-isoform-expression-rsem-tpm.polya.rds` | ||
|`pbta-isoform-expression-rsem-tpm.stranded.rds` || |
|`pbta-snv-consensus-mutation-tmb.tsv` | ||
|`pbta-snv-consensus-mutation.maf.tsv.gz` || |
|`pbta-snv-lancet.vep.maf.gz` | ||
|`pbta-snv-mutect2.vep.maf.gz` | ||
|`pbta-snv-strelka2.vep.maf.gz` | ||
|`pbta-snv-vardict.vep.maf.gz` | ||
|`pbta-sv-manta.tsv.gz`| ||
