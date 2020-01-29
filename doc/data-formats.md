# Data Formats in Data Download

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
A table with brief descriptions for each data file is provided in the `data-description.md` file included in the download.

## PBTA Data Files

PBTA data files are all files derived from samples (e.g., tumors, cell lines) that are processed upstream of this repository and are not the product of any analysis code in the `AlexsLemonade/OpenPBTA-analysis` repository.

### Quality Control Data

MendQC [output files](https://github.com/UCSC-Treehouse/mend_qc#output) `*readDist.txt` and `*bam_umend_qc.tsv`, along with a manifest mapping filename to biospecimen, are provided. 

 * `pbta-mend-qc-results.tar.gz`
 * `pbta-mend-qc-manifest.tsv`

### Somatic Single Nucleotide Variant (SNV) Data

Somatic Single Nucleotide Variant (SNV) data are provided in [Annotated MAF format](format/vep-maf.md) files for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-single-nucleotide-variant-calling) and denoted with the `pbta-snv` prefix.
A subset of TCGA brain tumor MAF files are also provided, along with a manifest, `pbta-tcga-manifest.tsv`, and denoted with the `pbta-tcga-snv` prefix.
Briefly, VCFs are VEP annotated and converted to MAF format via [vcf2maf](https://github.com/mskcc/vcf2maf/blob/master/maf2vcf.pl) to produce the following files:

* `pbta-snv-lancet.vep.maf.gz`
* `pbta-snv-mutect2.vep.maf.gz`
* `pbta-snv-strelka2.vep.maf.gz`
* `pbta-snv-vardict.vep.maf.gz`
* `pbta-tcga-snv-lancet.vep.maf.gz`
* `pbta-tcga-snv-mutect2.vep.maf.gz`
* `pbta-tcga-snv-strelka2.vep.maf.gz`

### Somatic Copy Number Variant (CNV) Data

Somatic Copy Number Variant (CNV) data are provided in a modified [SEG format](https://software.broadinstitute.org/software/igv/SEG) for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-copy-number-variant-calling) and denoted with the `pbta-cnv` prefix.
**Somatic copy number data is only generated for whole genome sequencing (WGS) samples.**

  * `pbta-cnv-cnvkit.seg.gz` is the the CNVkit SEG file. This file contains an additional column `copy.num` to denote copy number of each segment, derived from the CNS file output of the algorithm [described here](https://cnvkit.readthedocs.io/en/stable/fileformats.html).
  * `pbta-cnv-controlfreec.tsv.gz` is the ControlFreeC TSV file. It is a merge of `*_CNVs` files produced from the algorithm, and columns are [described here](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT).

#### A Note on Ploidy

The _copy number_ annotated in the CNVkit SEG file is annotated with respect to ploidy 2, however, the _status_ annotated in the ControlFreeC TSV file is annotated with respect to inferred ploidy from the algorithm, which is recorded in the `pbta-histologies.tsv` file. See the table below for examples of possible interpretations.
<br>
<br>

| Ploidy | Copy Number | Gain/Loss Interpretation     |
|--------|-------------|------------------------------|
| 2      | 0           | Loss; homozygous deletion    |
| 2      | 1           | Loss; hemizygous deletion    |
| 2      | 2           | Copy neutral                 |
| 2      | 3           | Gain; one copy gain          |
| 2      | 4           | Gain; two copy gain          |
| 2      | 5+          | Gain; possible amplification |
| 3      | 0           | Loss; 3 copy loss            |
| 3      | 1           | Loss; 2 copy loss            |
| 3      | 2           | Loss; 1 copy loss            |
| 3      | 3           | Copy neutral                 |
| 3      | 4           | Gain; one copy gain          |
| 3      | 5           | Gain; two copy gain          |
| 3      | 6+          | Gain; possible amplification |
| 4      | 0           | Loss; 4 copy loss            |
| 4      | 1           | Loss; 3 copy loss            |
| 4      | 2           | Loss; 2 copy loss            |
| 4      | 3           | Loss; 1 copy loss            |
| 4      | 4           | Copy neutral                 |
| 4      | 5           | Gain; one copy gain          |
| 4      | 6           | Gain; two copy gain          |
| 4      | 7+          | Gain; possible amplification |


### Somatic Structural Variant (SV) Data

Somatic Structural Variant Data (Somatic SV) are provided in the [Annotated Manta TSV](doc/format/manta-tsv-header.md) format produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-structural-variant-calling) and denoted with the `pbta-sv` prefix.
**Somatic structural variant data is only generated for whole genome sequencing (WGS) samples.**

* `pbta-sv-manta.tsv.gz`

### Gene Expression Data

Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a feature (e.g., gene or transcript) by sample matrix.
For each method/measure (e.g., RSEM TPM, RSEM isoform counts), there are two matrices provided: one for each library selection method (`poly-A`, `stranded`). 
See [this notebook](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/selection-strategy-comparison/01-selection-strategies.nb.html) for more information about _why_ this was necessary.
Gene expression are available in multiple forms in the following files:

* `pbta-gene-counts-rsem-expected_count.polya.rds`
* `pbta-gene-counts-rsem-expected_count.stranded.rds`
* `pbta-gene-expression-kallisto.polya.rds`
* `pbta-gene-expression-kallisto.stranded.rds`
* `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds`
* `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds`
* `pbta-gene-expression-rsem-fpkm-collapsed_table.polya.rds`
* `pbta-gene-expression-rsem-fpkm-collapsed_table.stranded.rds`
* `pbta-gene-expression-rsem-fpkm.polya.rds`
* `pbta-gene-expression-rsem-fpkm.stranded.rds`
* `pbta-gene-expression-rsem-tpm.polya.rds`
* `pbta-gene-expression-rsem-tpm.stranded.rds`
* `pbta-isoform-counts-rsem-expected_count.polya.rds`
* `pbta-isoform-counts-rsem-expected_count.stranded.rds`
* `pbta-isoform-expression-rsem-tpm.polya.rds`
* `pbta-isoform-expression-rsem-tpm.stranded.rds`

See [the data description file](data-description.md) for more information about the individual gene expression files.

If your analysis requires de-duplicated gene symbols as row names, please use the collapsed matrices provided as part of the data download ([see below](#collapsed-expression-matrices)).

STAR `*Log.final.out` files, along with a manifest mapping filename to biospecimen, are provided ([see section 4.1 here](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)). 

 * `pbta-star-log-final.tar.gz`
 * `pbta-star-log-manifest.tsv`

### Gene Fusion Data

Gene Fusions produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#rna-fusion-calling-and-prioritization) are provided as [Arriba TSV](format/arriba-tsv-header.md) and [STARFusion TSV](./format/starfusion-tsv-header.md) respectively.
These files are denoted with the prefix `pbta-fusion`.

* `pbta-fusion-arriba.tsv.gz`
* `pbta-fusion-starfusion.tsv.gz`


### Harmonized Clinical Data

[Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values in the following file:

* `pbta-histologies.tsv`

#### Mapping Between DNA-seq and RNA-seq Data for the Same Sample

Many analyses rely on examining DNA-seq (e.g., WGS/WXS/Panel) and RNA-seq data together.
The identifiers in the `Kids_First_Biospecimen_ID` column of the `pbta-histologies.tsv` file will be non-overlapping for different experimental strategies, as these identifiers map to a _library or assay_.
**Analysts should use the identifiers in the `sample_id` to connect DNA-seq and RNA-seq assays from the same sample.**

For an example, see [the `sample_id` mapping step of the OncoPrint pipeline](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/oncoprint-landscape/00-map-to-sample_id.R).

Note that some individual participants (tracked via the `Kids_First_Participant_ID` column) will have multiple samples included in the PBTA dataset.
Please [see the independent specimens files section](#independent-specimen-lists) for more information.

<!-- TODO: add a notebook with an example for how to map between identifiers -->

## Analysis Files

Analysis files are created by a script in `analyses/*`. 
They can be viewed as _derivatives_ of PBTA data files.

### Consensus Mutation Files

Consensus mutation files are products of the [`analyses/snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) analysis module.

  * `pbta-snv-consensus-mutation.maf.tsv.gz` contains consensus calls for SNVs and small indels. The calls in the file are created as the intersection of calls from Strelka2, Mutect2, Lancet, where position, change, and sample were the same for all callers.
  Multinucleotide variant calls from Mutect2 and Lancet were separated into consecutive SNVs before merging.
  All columns in the included file are derived from the Strelka2 calls.
  Note that this file is not strictly a MAF file, as it adds a Variant Allele Frequency (`VAF`) column and does not contain a version comment as the first line.
  * `pbta-snv-consensus-mutation-tmb-all.tsv` includes tumor mutation burden statistics that are calculated based calculated from Strelka2 and Mutect2 SNV consensus, and the intersection of Strelka2 and Mutect2 BED windows sizes.
  * `pbta-snv-consensus-mutation-tmb-coding.tsv` contains coding only tumor mutation burden statistics calculated from the number of coding sequence Strelka2, Mutect2, and Lancet consensus SNVs and size of the intersection of all three callers' BED windows and the Gencode v27 coding sequence regions. 

### Collapsed Expression Matrices

Collapsed expression matrices are products of the [`analyses/collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) analysis module.
In cases where more than one Ensembl gene identifier maps to the same gene symbol, the instance of the gene symbol with the maximum mean FPKM in the RSEM FPKM file is retained to produce the following files:

* `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds`  
* `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds`

### Independent Specimen Lists

Independent specimen lists are the products of the [`analyses/independent-samples`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/independent-samples) analysis module.
For participants with multiple tumor specimens, independent specimen lists are provided as TSV files with columns for participant ID (`Kids_First_Participant_ID`) and specimen ID (`Kids_First_Biospecimen_ID`).
The methods are [described here](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples)). 
These files are used for analyses such as mutation co-occurence, where repeated samples might cause bias. 
  
  * `independent-specimens.wgs.primary.tsv` with WGS samples and only primary tumors.
  * `independent-specimens.wgs.primary-plus.tsv` as above, but including non-primary tumors where a primary tumor sample is not available.
  * `independent-specimens.wgswxs.primary.tsv` Only primary tumors, but with WXS where WGS is not available. 
  * `independent-specimens.wgswxs.primary-plus.tsv` as above, but including non-primary tumors where a primary tumor sample is not available.

Note that these independent specimen files do not address the issue of participants with multiple tumor specimens in RNA-seq data at this time.

### Derived Fusion Files

The filtered and prioritized fusion and downstream files are a product of the [`analyses/fusion_filtering`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering) analysis module. 

  * `pbta-fusion-putative-oncogenic.tsv` contains the filtered and prioritized fusions. 
  The methods are [described here](https://alexslemonade.github.io/OpenPBTA-manuscript/#fusion-prioritization).
  * `pbta-fusion-recurrently-fused-genes-byhistology.tsv` is a table that includes counts of recurrently fused genes by broad histology.
  * `pbta-fusion-recurrently-fused-genes-bysample.tsv` contains a binary matrix that denotes the presence or absence of a recurrently fused gene in an individual RNA-seq specimen.
  * `fusion_summary_embryonal_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent embryonal tumor fusions of interest per individual RNA-seq specimen.
  * `fusion_summary_ependymoma_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent ependymal tumor fusions of interest per individual RNA-seq specimen.

### Copy Number Files

`pbta-cnv-cnvkit-gistic.zip` is the output of running GISTIC 2.0 on the CNVkit results (`pbta-cnv-cnvkit.seg`).
The script used to run GISTIC can be [found here](https://github.com/d3b-center/publication_workflows/blob/master/openPBTA/run-gistic.sh).

### Consensus Copy Number File

Copy number consensus calls from the copy number and structural variant callers are a product of the [`analyses/copy_number_consensus_call`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/copy_number_consensus_call) analysis module. 
 * `cnv_consensus.tsv` contains consensus regions from two or more callers, with columns described in the [analysis README](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/copy_number_consensus_call/README.md).

## Data Caveats

The clinical manifest will be updated and versioned as molecular subgroups are identified based on genomic analyses.
