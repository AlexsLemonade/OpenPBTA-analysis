# Data Formats in Data Download

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
A table with brief descriptions for each data file is provided in the `data-files-description.md` file included in the download.

## Processed Data Files

Processed data files are all files derived from samples (e.g., tumors, cell lines) that are processed upstream of this repository and are not the product of any analysis code in the `AlexsLemonade/OpenPBTA-analysis` repository.


### Somatic Single Nucleotide Variant (SNV) Data

Somatic Single Nucleotide Variant (SNV) data are provided in [Annotated MAF format](format/vep-maf.md) files for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-single-nucleotide-variant-calling) and denoted with the `pbta-snv` prefix. Briefly, VCFs are VEP annotated and converted to MAF format via [vcf2maf](https://github.com/mskcc/vcf2maf/blob/master/maf2vcf.pl) to produce the following files.

* `kfnbl-snv-lancet.vep.maf.gz`
* `kfnbl-snv-mutect2.vep.maf.gz`
* `kfnbl-snv-strelka2.vep.maf.gz`
* `kfnbl-snv-vardict.vep.maf.gz`

### Consensus Somatic Variant Data
Somatic calls that are retained if they are supported by atleast 2 callers OR marked as `HotSpotAllele` because they overlap SNV/INDELs considered as [Cancer Hotspots](https://www.cancerhotspots.org/#/download) OR are TERT promoter SNVs. Please find additional information [here](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md) 

* `snv-consensus-plus-hotspots.maf.tsv.gz`

### Gene Expression Data

Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a feature (e.g., gene or transcript) by sample matrix.
For each method/measure (e.g., RSEM TPM, RSEM isoform counts), there are two matrices provided: one for each library selection method (`poly-A`, `stranded`). 
See [this notebook](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/selection-strategy-comparison/01-selection-strategies.nb.html) for more information about _why_ this was necessary.
Gene expression are available in multiple forms in the following files:

* `gene-counts-rsem-expected_count.rds`
* `gene-expression-rsem-tpm.rds`

See [the data description file](data-description.md) for more information about the individual gene expression files.

If your analysis requires de-duplicated gene symbols as row names, please use the collapsed matrices provided as part of the data download ([see below](#collapsed-expression-matrices)).

### Gene Fusion Data

Gene Fusions produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#rna-fusion-calling-and-prioritization) are provided as [Arriba TSV](format/arriba-tsv-header.md) and [STARFusion TSV](./format/starfusion-tsv-header.md) respectively.
These files are denoted with the prefix `pbta-fusion`.

* `fusion-arriba.tsv.gz`
* `fusion-starfusion.tsv.gz`


### Harmonized Clinical Data

[Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values in the following file:

* `histologies.tsv`

## Analysis Files

Analysis files are created by a script in `analyses/*`. 
They can be viewed as _derivatives_ of Processed data files.

### Collapsed Expression Matrices

Collapsed expression matrices are products of the [`analyses/collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) analysis module.
In cases where more than one Ensembl gene identifier maps to the same gene symbol, the instance of the gene symbol with the maximum mean FPKM in the RSEM FPKM file is retained to produce the following files:

* `gene-counts-rsem-expected_count-collapsed.rds`  
* `gene-counts-rsem-expected_count-collapsed.rds`
