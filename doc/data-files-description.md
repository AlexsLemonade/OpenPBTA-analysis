## Data file descriptions

This document contains information about all data files associated with this project. Each file will have the following association information:

+ **File type** will be one of:
	+ *Reference file*: Obtained from an external source/database. When known, the obtained data and a link to the external source is included.
	+ *Modified reference file*: Obtained from an external source/database but modified for OpenPBTA use. 
	+ *Processed data file*: Data that are processed upstream of the analysis project, e.g., the output of a somatic single nucleotide variant method. Links to the relevant D3B Center or Kids First workflow (and version where applicable) are included in **Origin**.
	+ *Analysis file*: Any file created by a script in `analyses/*`. 
+ **Origin**
	+ For _Processed data files_, a link the relevant D3B Center or Kids First workflow (and version where applicable).
	+ When applicable, a link to the specific *script* that produced (or modified, for *Modified reference file* types) the data.
+ **File description**
	+ A *brief* one sentence description of what the file contains (e.g., bed files contain coordinates for features XYZ).

### current release (v4)

| **File name** |  **File Type** | **Origin** | **File Description** |
|---------------|----------------|------------------------|-----------------------|
|`intersect_cds_lancet_strelka_mutect_WGS.bed` | Analysis file | [`analyses/snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/) | Intersection of `gencode.v27.primary_assembly.annotation.gtf.gz` CDS with Lancet, Strelka2, Mutect2 regions
|`intersect_strelka_mutect_WGS.bed` | Analysis file | [`analyses/snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/) | Intersection of `gencode.v27.primary_assembly.annotation.gtf.gz` CDS with Strelka2 and Mutect2 regions called
|`fusion-arriba.tsv.gz` | Processed data file | [Gene fusion detection](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#gene-fusion-detection); [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_RNAseq_workflow.cwl) | Fusion - [Arriba TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/arriba-tsv-header.md), annotated with FusionAnnotator
|`fusion-starfusion.tsv.gz` | Processed data file | [Gene fusion detection](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#gene-fusion-detection); [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_RNAseq_workflow.cwl) | Fusion - [STARFusion TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/starfusion-tsv-header.md)
|`gene-counts-rsem-expected_count-collapsed.rds` | Analysis file | PBTA+GMKF[`analysis/collapse-rnaseq`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq) ;GTEx v8 release| Gene expression - RSEM expected_count for each samples collapsed to gene symbol (gene-level)
|`gene-expression-rsem-tpm-collapsed.rds` | Analysis file | PBTA+GMKF[`analysis/collapse-rnaseq`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq);GTEx v8 release | Gene expression - RSEM TPM for each samples collapsed to gene symbol (gene-level)
|`kfnbl-snv-lancet.vep.maf.gz` | Processed data file | [Somatic mutation calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-mutation-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc-lancet-wf.cwl) | Somatic SNV - Lancet [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`kfnbl-snv-mutect2.vep.maf.gz` | Processed data file | [Somatic mutation calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-mutation-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc_strelka2_mutect2_manta_workflow.cwl) | Somatic SNV - Mutect2 [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`kfnbl-snv-strelka2.vep.maf.gz` | Processed data file | [Somatic mutation calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-mutation-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc_strelka2_mutect2_manta_workflow.cwl) | Somatic SNV - Strelka2 [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`kfnbl-snv-vardict.vep.maf.gz` | Processed data file | [Somatic mutation calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-mutation-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc-vardict-wf.cwl) | Somatic SNV - VarDict [annotated MAF file](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/vep-maf.md)
|`WGS.hg38.lancet.300bp_padded.bed` | Reference Target/Baits File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | WGS.hg38.lancet.unpadded.bed file with each region padded by 300 bp
|`WGS.hg38.lancet.unpadded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) |  hg38 WGS regions created using UTR, exome, and start/stop codon features of the GENCODE 31 reference, augmented with PASS variant calls from Strelka2 and Mutect2
|`WGS.hg38.mutect2.vardict.unpadded.bed` | Reference Regions File  | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) |  hg38 BROAD Institute interval calling list (restricted to Chr1-22,X,Y,M and non-N regions) used for Mutect2 and VarDict variant callers
|`WGS.hg38.strelka2.unpadded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | hg38 BROAD Institute interval calling list (restricted to Chr1-22,X,Y,M) used for Strelka2 variant caller
|`WGS.hg38.vardict.100bp_padded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | `WGS.hg38.mutect2.vardict.unpadded.bed` with each region padded by 100 bp used for VarDict variant caller
| `snv-consensus-plus-hotspots.maf.tsv.gz` | Processed data file | [Consensus calling method](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md) | Consensus (2 of 4) maf for PBTA + GMKF) |
