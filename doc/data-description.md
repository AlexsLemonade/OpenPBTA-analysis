## Data file descriptions

This document contains information about all data files associated with this project. Each file will have the following association information:

+ **File type** will be one of:
	+ *Reference file*: Obtained from an external source/database. When known, the obtained data and a link to the external source is included.
	+ *Modified reference file*: Obtained from an external source/database but modified for OpenPBTA use. 
	+ *PBTA data file*: Pediatric Brain Tumor Atlas data that are processed upstream of the OpenPBTA project, e.g., the output of a somatic single nucleotide variant method. Links to the relevant D3B Center or Kids First workflow (and version where applicable) are included in **Origin**.
	+ *Analysis file*: Any file created by a script in `analyses/*`. 
+ **Origin**
	+ For _PBTA data files_, a link the relevant D3B Center or Kids First workflow (and version where applicable).
	+ When applicable, a link to the specific *script* that produced (or modified, for *Modified reference file* types) the data.
+ **File description**
	+ A *brief* one sentence description of what the file contains (e.g., bed files contain coordinates for features XYZ).



### current release (release-v11-20191126)

| **File name** |  **File Type** | **Origin** | **File Description** |
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
