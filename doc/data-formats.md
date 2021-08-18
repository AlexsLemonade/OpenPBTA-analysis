# Data Formats in Data Download

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
A table with brief descriptions for each data file is provided in the `data-files-description.md` file included in the download.

## Processed Data Files

Processed data files are all files derived from samples (e.g., tumors, cell lines) that are processed upstream of this repository and are not the product of any analysis code in the `AlexsLemonade/OpenPBTA-analysis` repository.


### Somatic Single Nucleotide Variant (SNV) Data

Somatic Single Nucleotide Variant (SNV) data are provided in [Annotated MAF format](format/vep-maf.md) files for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-single-nucleotide-variant-calling) and denoted with the `snv` prefix. Briefly, VCFs are VEP annotated and converted to MAF format via [vcf2maf](https://github.com/mskcc/vcf2maf/blob/master/maf2vcf.pl) to produce the MAF files. 

* `snv-lancet.vep.maf.gz`
* `snv-mutect2.vep.maf.gz`
* `snv-strelka2.vep.maf.gz`
* `snv-vardict.vep.maf.gz`

### Consensus Somatic Variant Data
Somatic calls that are retained if they are supported by atleast 2 callers OR marked as `HotSpotAllele` because they overlap SNV/INDELs considered as [Cancer Hotspots](https://www.cancerhotspots.org/#/download) OR are TERT promoter SNVs. Please find additional information [here](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md) 

* `snv-consensus-plus-hotspots.maf.tsv.gz`

### Somatic Copy Number Variant (CNV) Data

Somatic Copy Number Variant (CNV) data are provided in a modified [SEG format](https://software.broadinstitute.org/software/igv/SEG) for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-copy-number-variant-calling) and denoted with the `cnv` prefix.
**Somatic copy number data is only generated for whole genome sequencing (WGS) samples.**

  * `cnv-cnvkit.seg.gz` is the the CNVkit SEG file. This file contains an additional column `copy.num` to denote copy number of each segment, derived from the CNS file output of the algorithm [described here](https://cnvkit.readthedocs.io/en/stable/fileformats.html).
  * `cnv-controlfreec.tsv.gz` is the ControlFreeC TSV file. It is a merge of `*_CNVs` files produced from the algorithm, and columns are [described here](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT).

#### A Note on Ploidy

The _copy number_ annotated in the CNVkit SEG file is annotated with respect to ploidy 2, however, the _status_ annotated in the ControlFreeC TSV file is annotated with respect to inferred ploidy from the algorithm, which is recorded in the `histologies.tsv` file. 

### Gene Expression Data

Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a feature (e.g., gene or transcript) by sample matrix.
Gene expression are available in multiple forms in the following files:

* `gene-counts-rsem-expected_count.rds`
* `gene-expression-rsem-tpm.rds`

See [the data description file](data-description.md) for more information about the individual gene expression files.

If your analysis requires de-duplicated gene symbols as row names, please use the collapsed matrices provided as part of the data download ([see below](#collapsed-expression-matrices)).

### Gene Fusion Data

Gene Fusions produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#rna-fusion-calling-and-prioritization) are provided as [Arriba TSV](format/arriba-tsv-header.md) and [STARFusion TSV](./format/starfusion-tsv-header.md) respectively.
These files are denoted with the prefix `fusion`.

* `fusion-arriba.tsv.gz`
* `fusion-starfusion.tsv.gz`

### Structural Variant Data

Structural Variants data produced by the [`MANTA` package](https://pubmed.ncbi.nlm.nih.gov/26647377/; https://github.com/Illumina/manta) is 
provided as:

* `sv-manta.tsv.gz`

### Harmonized Clinical Data

[Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values in the following file:

* `histologies.tsv`

### Independent Sample List

[Independent sample list](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/independent-samples) are released as tab separated values in the following file. `wgswxspanel` indicates it includes all experimental strategy for DNA sequencing, `RNA` includes all RNA samples and `primary` and `relapse` specifies which tumor type the samples are from. `eachcohort` indicates the selection is cohort-based.

* `independent-specimens.wgswxspanel.primary.eachcohort.tsv`
* `independent-specimens.wgswxspanel.relapse.eachcohort.tsv`
* `independent-specimens.rnaseq.primary.eachcohort.tsv`
* `independent-specimens.rnaseq.relapse.eachcohort.tsv`

## Analysis Files

Analysis files are created by a script in `analyses/*`. 
They can be viewed as _derivatives_ of Processed data files.

### Collapsed Expression Matrices

Collapsed expression matrices are products of the [`analyses/collapse-rnaseq`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq) analysis module.
In cases where more than one Ensembl gene identifier maps to the same gene symbol, the instance of the gene symbol with the maximum mean FPKM in the RSEM FPKM file is retained to produce the following files:

* `gene-counts-rsem-expected_count-collapsed.rds`  
* `gene-expression-rsem-tpm-collapsed.rds`

### Derived Copy Number Files

#### Consensus Copy Number File

Copy number consensus calls from the copy number and structural variant callers are a product of the [`analyses/copy_number_consensus_call`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call) analysis module. 

* `cnv-consensus.seg.gz` contains consensus segments and segment means (log R ratios) from two or more callers, as described in the [analysis README](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call/README.md) - contains only WGS samples. 

##### Focal Copy Number Files

Focal copy number files map the consensus calls (genomic segments) in WGS samples to genes for downstream analysis and are a product of the [`analysis/focal-cn-file-preparation`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation).
Note: these files contain biospecimens and genes with copy number changes. 

  - `consensus_seg_annotated_cn_autosomes.tsv.gz` contains focal gene copy number alterations for all autosomes.
  - `consensus_seg_annotated_cn_x_and_y.tsv.gz` contains focal gene copy number alterations for the sex chromosomes.
  
  Focal copy number files in WXS samples only uses results from CNVkit and no consensus calling is required. [`analysis/focal-cn-file-preparation`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation).
Note: these files contain biospecimens and genes with copy number changes. 

  - `cnvkit_annotated_cn_wxs_autosomes.tsv.gz` contains focal gene copy number alterations for all autosomes.
  - `cnvkit_annotated_cn_wxs_x_and_y.tsv.gz` contains focal gene copy number alterations for the sex chromosomes.
  
  Additionally, autosomes file and x_and_y file for either WGS or WXS to generate the two combined files as followed:
  - `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz`
  - `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz` 
  And these two files are further merged to generate:
  - `consensus_wgs_plus_cnvkit_wxs.tsv.gz`
