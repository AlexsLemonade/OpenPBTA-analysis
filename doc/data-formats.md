## Data Formats

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
A table with brief descriptions for each data file is provided in the `data-description.md` 

* Somatic Single Nucleotide Variant (SNV) data are provided in [Annotated MAF format](doc/format/vep-maf.md) files for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-single-nucleotide-variant-calling).
  * Consensus calls for SNVs and small indels in the file `pbta-snv-consensus-mutation.maf.tsv.gz` are created as the intersection of calls from Strelka2, Mutect2, Lancet, where position, change, and sample were the same for all callers.
  Multinucleotide variant calls from Mutect2 and Lancet were separated into consecutive SNVs before merging.
  All columns in the included file are derived from the Strelka2 calls.
  Note that this file is not strictly a MAF file, as it adds a Variant Allele Frequency (`VAF`) column and does not contain a version comment as the first line.
  * Tumor mutation burden statistics are calculated based on the mutations included in `pbta-snv-consensus-mutation.maf.tsv.gz` file and by using Strelka2 counts and BED window sizes.
  These values are saved to `pbta-snv-consensus-mutation-tmb.tsv`
* Somatic Copy Number Variant (CNV) data are provided in a modified [SEG format](https://software.broadinstitute.org/software/igv/SEG) for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-copy-number-variant-calling).
  * The CNVkit SEG file has an additional column `copy.num` to denote copy number of each segment, derived from the CNS file output of the algorithm described [here](https://cnvkit.readthedocs.io/en/stable/fileformats.html).
  * The ControlFreeC TSV file is a merge of `*_CNVs` files produced from the algorithm, and columns are described [here](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT).
  * NOTE: The _copy number_ annotated in the CNVkit SEG file is annotated with respect to ploidy 2, however, the _status_ annotated in the ControlFreeC TSV file is annotated with respect to inferred ploidy from the algorithm, which is recorded in the `pbta_histologies.tsv` file. See the table below for examples of possible interpretations.

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

* Somatic Structural Variant Data (Somatic SV) are provided in the [Annotated Manta TSV](doc/format/manta-tsv-header.md) format produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-structural-variant-calling).
* Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a gene by sample matrix.
	* If your analysis requires de-duplicated gene symbols as row names, please use the collapsed matrices provided as part of the data download (`pbta-gene-expression-rsem-fpkm-collapsed.polya.rds`, `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds`).
* Gene Fusions produced by the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#rna-fusion-calling-and-prioritization) are provided as [Arriba TSV](doc/format/arriba-tsv-header.md) and [STARFusion TSV](doc/format/starfusion-tsv-header.md) respectively.
* [Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values.
* For participants with multiple tumor specimens, [Independent specimen lists](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) are provided as TSV files with columns for participant ID and specimen ID. 
These files are used for analyses such as mutation co-occurence, where repeated samples might cause bias.
There are four of these files:
  1. `independent-specimens.wgs.primary.tsv` with WGS samples and only primary tumors
  2. `independent-specimens.wgs.primary-plus.tsv` as above, but including non-primary tumors where a primary tumor sample is not available
  3. `independent-specimens.wgswxs.primary.tsv` Only primary tumors, but with WXS where WGS is not available 
  4. `independent-specimens.wgswxs.primary-plus.tsv` as above, but including non-primary tumors where a primary tumor sample is not available.
  

### Data Caveats

The clinical manifest will be updated and versioned as molecular subgroups are identified based on genomic analyses.
