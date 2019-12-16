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
|`GRCh38.primary_assembly.genome.fa.gz` | Reference file | GENCODE v27 | hg38 primary assembly genome sequence FASTA file
|`StrexomeLite_Targets_CrossMap_hg38_filtered_chr_prefixed.bed` | | |
|`StrexomeLite_hg38_liftover_100bp_padded.bed`| | |
|`WGS.hg38.lancet.300bp_padded.bed` | | |
|`WGS.hg38.lancet.unpadded.bed` | | |
|`WGS.hg38.mutect2.unpadded.bed` | | |
|`WGS.hg38.strelka2.unpadded.bed` | ||
|`WGS.hg38.vardict.100bp_padded.bed` | ||
|`WXS.hg38.100bp_padded.bed` | ||
|`gencode.v27.primary_assembly.annotation.gtf.gz` | Reference file | GENCODE v27 | hg38 gene annotation on primary assembly (reference chromosomes and scaffolds)
|`independent-specimens.wgs.primary-plus.tsv` | Analysis file |[`analyses/independent-samples`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/independent-samples)| Independent specimens list for WGS sample, primary + non-primary when no primary sample is available
|`independent-specimens.wgs.primary.tsv` | Analysis file | [`analyses/independent-samples`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/independent-samples) | Independent specimens list for WGS samples, primary only
|`independent-specimens.wgswxs.primary-plus.tsv` | Analysis file | [`analyses/independent-samples`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/independent-samples) | Independent specimens list for WGS and WXS samples, primary + non-primary when no primary sample is available
|`independent-specimens.wgswxs.primary.tsv` | Analysis file | [`analyses/independent-samples`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/independent-samples) | Independent specimens list for WGS and WXS samples, primary only
|`pbta-cnv-cnvkit.seg.gz` | PBTA data file || Somatic Copy Number Variant - CNVkit [SEG file](https://cnvkit.readthedocs.io/en/stable/fileformats.html#seg)
|`pbta-cnv-controlfreec.tsv.gz` | PBTA data file || Somatic Copy Number Variant - TSV file that is a merge of [ControlFreeC `*_CNVs` files](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT)
|`pbta-fusion-arriba.tsv.gz` | PBTA data file || Fusion - [Arriba TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/arriba-tsv-header.md)
|`pbta-fusion-putative-oncogenic.tsv` | Analysis file | [`analyses/fusion_filtering`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering) | Filtered and prioritized fusions 
|`pbta-fusion-starfusion.tsv.gz` | PBTA data file || Fusion - [STARFusion TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/starfusion-tsv-header.md)
|`pbta-gene-counts-rsem-expected_count.polya.rds` | PBTA data file || Gene expression - RSEM expected counts for poly-A samples (gene-level)
|`pbta-gene-counts-rsem-expected_count.stranded.rds` | PBTA data file | | Gene expression - RSEM  expected counts for stranded samples (gene-level)
|`pbta-gene-expression-kallisto.polya.rds` | PBTA data file || Gene expression - kallisto TPM for poly-A samples (transcript-level)
|`pbta-gene-expression-kallisto.stranded.rds` | PBTA data file | | Gene expression - kallisto TPM for stranded samples (transcript-level)
|`pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` | Analysis file | [`analyses/collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) | Gene expression - RSEM FPKM for poly-A samples collapsed to gene symbol (gene-level)
|`pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` | Analysis file | [`analyses/collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) | Gene expression - RSEM FPKM for stranded samples collapsed to gene symbol (gene-level)
|`pbta-gene-expression-rsem-fpkm.polya.rds` | PBTA data file || Gene expression - RSEM FPKM for poly-A samples (gene-level)
|`pbta-gene-expression-rsem-fpkm.stranded.rds` | PBTA data file | | Gene expression - RSEM FPKM for stranded samples (gene-level)
|`pbta-gene-expression-rsem-tpm.polya.rds` | PBTA data file || Gene expression - RSEM TPM for poly-A samples (gene-level)
|`pbta-gene-expression-rsem-tpm.stranded.rds` | PBTA data file | | Gene expression -RSEM TPM for stranded samples (gene-level)
|`pbta-histologies.tsv` | PBTA data file || Harmonized clinical metadata file (see data dictionary [here](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#clinical-data-harmonization))
|`pbta-isoform-counts-rsem-expected_count.polya.rds` | PBTA data file || Gene expression -RSEM expected counts for poly-A samples (transcript-level)
|`pbta-isoform-counts-rsem-expected_count.stranded.rds` | PBTA data file | |Gene expression - RSEM expected counts for stranded samples (transcript-level)
|`pbta-isoform-expression-rsem-tpm.polya.rds` | PBTA data file || Gene expression - RSEM TPM for poly-A samples (transcript-level)
|`pbta-isoform-expression-rsem-tpm.stranded.rds` | PBTA data file | | Gene expression - RSEM TPM for stranded samples (transcript-level)
|`pbta-snv-consensus-mutation-tmb.tsv` | Analysis file | [`analyses/snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/) | Tumor mutation burden statistics calculated from consensus SNV, using Strelka2 counts and BED window sizes
|`pbta-snv-consensus-mutation.maf.tsv.gz` | Analysis file | [`analyses/snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/)  | Consensus calls for SNVs and small indels; columns in the included file are derived from the Strelka2.
|`pbta-snv-lancet.vep.maf.gz` | PBTA data file | | Somatic SNV - Lancet [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`pbta-snv-mutect2.vep.maf.gz` | PBTA data file || Somatic SNV - Mutect2 [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`pbta-snv-strelka2.vep.maf.gz` | PBTA data file || Somatic SNV - Strelka2 [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`pbta-snv-vardict.vep.maf.gz` | PBTA data file || Somatic SNV - VarDict [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`pbta-sv-manta.tsv.gz`| PBTA data file || Somatic Structural Variant - Manta output, annotated with AnnotSV
