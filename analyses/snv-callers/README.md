# SNV caller comparison analysis

This analysis evaluates [MAF files](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)
from different SNV callers and compares their output.
The GDC has [good documentation on the fields](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) contained in a standard MAF file.

**Table of Contents**
* [Individual caller evaluation](#individual-caller-evaluation)  
  * [Base change analysis](#base-change-analysis)  
  * [Variant allele fraction calculation](#variant-allele-fraction-calculation)  
  * [Genomic regional analyses](#genomic-regional-analyses)  
  * [Tumor mutation burden](#tumor-mutation-burden-calculation)
  * [COSMIC mutation overlap](#cosmic-mutation-overlap)  
* Comparison analysis - in development
 * [Mutation IDs](#mutation-ids)
* [Overall file structure](#overall-file-structure)
* [Summary of functions](#summary-of-special-functions)
## Individual Caller Evaluation

The first step in this analysis is an individual evaluation of each MAF from
each caller.
This analysis prints out a report that is further subdivided by the WGS and WXS
samples.

Overall reports for each caller and strategy can be found:
* `results/<caller_name>/<caller_name>_wgs_report.html`
* `results/<caller_name>/<caller_name>_ws_report.html`

### Base change analysis

To evaluate base changes I summarized standard MAF fields as two new variables:
The `base_change` variable that indicates the exact change in bases from
concatenating `Reference_Allele`, `>`, and `Allele`.
The `change` variable is made from the `base_change` variable but groups
together deletions, insertions, and long (more than a SNV) as their own groups.

*Output for this analysis*
* `results/<caller_name>/<caller_name>_vaf.tsv`
* `plots/<caller_name>/<caller_name>_<strategy>_base_change.png`

### Variant Allele Fraction Calculation

Calculate variant allele fraction (VAF) for each variant.
```
vaf = (t_alt_count) / (t_ref_count + t_alt_count)
```
 This is following the [code used in
`maftools`](https://github.com/PoisonAlien/maftools/blob/1d0270e35c2e0f49309eba08b62343ac0db10560/R/plot_vaf.R#L39).
The VAF calculations and other special variables are added to the MAF fields
 and written to a TSV ending in `_vaf.tsv` in the caller's results folder.

 *Output for this analysiss*
 * `results/<caller_name>/<caller_name>_vaf.tsv`
 * `plots/<caller_name>/<caller_name>_<strategy>_depth_vs_vaf.png`

### Genomic Regional Analyses

To analyze what genomic regions the variants are from, I used [Annotatr
package](https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html) to obtain hg38 genome annotations. This Annotatr object is stored as an RDS file: `hg38_genomic_region_annotations.rds` in the `scratch` directory.
Mutations are assigned all annotations that they overlap (using  
  `GenomicRanges::overlap`).

*Output for this analysis*
* `results/<caller_name>/<caller_name>_regions.tsv`
* `plots/<caller_name>/<caller_name>_<strategy>_snv_regions.png`

### Tumor Mutation Burden Calculation

To calculate TMB, WXS or WGS BED regions are used as the denominator, depending on the sample's processing strategy.

```
TMBwxs = sum(mutation_w-in_bedwxs)/(wxs_genome_size/1000000)
TMBwgs = sum(all_mutations)/(wgs_genome_size/1000000)
```

Where genome size is calculated from the respective BED file as:
```
genome_size = sum(End_Position - Start_Position)
```

BED regions for WXS samples can be [found here](https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/Exome-AZ_V2.bed).
BED regions used for WGS samples are caller specific are from <unknown as of now>
The sample-wise TMB calculations written to a TSV ending in `_tmb.tsv` in the
caller's results folder.

*Output for this analysis*
* `results/<caller_name>/<caller_name>_tmb.tsv`
* `plots/<caller_name>/<caller_name>_<strategy>_tmb_plot.png`

### COSMIC Mutation Overlap

The COSMIC mutation data were obtained from https://cancer.sanger.ac.uk/cosmic/download
*To run this analysis, you need to obtain these data.*
The full, unfiltered somatic mutations file `CosmicMutantExport.tsv` for grch38
is used here and the genomic coordinates is arranged to be in BED format.
COSMIC mutations are overlapped with the present data's mutations using
GenomicRanges. The outcome of this overlap is added to the VAF data.frame with
two TRUE/FALS columns: `overlap_w_cosmic` is TRUE for mutations that overlap
with COSMIC mutations, while `same_as_cosmic` is TRUE when the base change
summary is also identical. The VAF for mutations that are or are not overlapping
with COSMIC mutations are then plotted in a violin plot.

*Output for this analysis*
* `results/<caller_name>/<caller_name>_vaf.tsv`
* `plots/<caller_name>/<caller_name>_<strategy>_cosmic_plot.png`

## Comparison of Callers - coming soon

### Mutation IDs  

In order to compare mutations across callers, I created a `mutation_id` from
combining information from standard MAF fields.

`mutation_id` is a concatenation of:  
* `Hugo_Symbol`  
* [`change`](#base-change-analysis)  
* `Start_Position`  
* `Tumor_Sample_Barcode` (the sample ID)  

If mutation_id's are identical among MAF files, they are considered the same.

## Overall file structure
```
OpenPBTA-analysis
├── analyses
│   └── snv-callers
│       ├── run_variant_caller_reports.R
│       ├── scripts
│       │   ├── 00-set-up.R
│       │   ├── 01-calculate_vaf_tmb.R
│       │   └── 02-run_eval.R
│       ├── functions
│       │    ├── plot_functions.R
│       │    └── wrangle_functions.R
│       ├── bed_regions
│       │   ├── wxs_bed_regions_all.tsv
│       │   ├── lancet_wgs_bed_regions.tsv
│       │   ├── mutect2_wgs_bed_regions.tsv
│       │   ├── strelka2_wgs_bed_regions.tsv
│       │   └── vardict_wgs_bed_regions.tsv
│       ├── results
│       │   ├── lancet
│       │   │   ├── lancet_vaf.tsv
│       │   │   ├── lancet_vaf.tsv.zip
│       │   │   ├── lancet_tmb.tsv
│       │   │   ├── lancet_tmb.tsv.zip
│       │   │   ├── lancet_region.tsv
│       │   │   ├── lancet_region.tsv.zip
│       │   │   ├── lancet_wxs_report.html
│       │   │   ├── lancet_wxs_report.Rmd
│       │   │   ├── lancet_wgs_report.html
│       │   │   ├── lancet_wgs_report.Rmd
│       │   │   └── lancet_metadata_filtered.tsv
│       │   ├── mutect2
│       │   │   └── ...
│       │   ├── strelka2
│       │   │   └── ...
│       │   └── vardict
│       │       └── ...
│       ├── plots
│       │   ├── lancet
│       │   │   ├── lancet_wgs_base_change.png
│       │   │   ├── lancet_wgs_cosmic_plot.png
│       │   │   ├── lancet_wgs_depth_vs_vaf.png
│       │   │   ├── lancet_wgs_snv_region.png
│       │   │   ├── lancet_wgs_tmb_plot.png
│       │   │   ├── lancet_wxs_base_change.png
│       │   │   ├── lancet_wxs_cosmic_plot.png
│       │   │   ├── lancet_wxs_depth_vs_vaf.png
│       │   │   ├── lancet_wxs_snv_region.png
│       │   │   └── lancet_wxs_tmb_plot.png
│       │   ├── mutect2
│       │   │   └── ...
│       │   ├── strelka2
│       │   │   └── ...
│       │   └── vardict
│       │       └── ...
│       └── template
│           └── variant_caller_report_template.Rmd
├── data
├── scratch
│    ├── cosmic_variants_cleaned.tsv
│    └── hg38_genomic_region_annotations.rds
```

## Summary of special functions

### Wrangling functions
|Function Name|Output created|Main Arguments|
|-------------|--------------|---------|
|`calculate_vaf`|Columns: VAF, mutation_id, base_change, change| A MAF formatted data.frame|
|`maf_to_granges`|A `GenomicRanges` format object|A MAF formatted data.frame|
|`wxs_bed_filter`|A filtered MAF df with only mutations within the provided BED regions|A MAF formatted data.frame and a BED formatted data.frame|
|`calculate_tmb`|A data.frame with TMB stat per sample |A `wxs_bed_filter`ed MAF formatted data.frame, a WGS and WXS genome sizes|
|`annotr_maf`|A data.frame with genomic annotations for the provided MAF data.frame|A MAF formatted data.frame and a [built AnnotatR annotation object](https://rdrr.io/bioc/annotatr/man/build_annotations.html)|
|`find_cosmic_overlap`|Find overlap with [COSMIC](https://cancer.sanger.ac.uk/cosmic) mutations|A MAF formatted data.frame with `change` and variable from the `calculate_vaf` function|

### Plotting functions

All plotting functions use `ggplot2` and all use a argument: `exp_strategy` that
determines whether to plot WGS or WXS samples only or both. Default is to plot
both.

|Function Name|Plot output|Main Arguments|
|-------------|-----------|---------|
|`base_change_plot`|Base change ggplot barplot|A MAF formatted data.frame with the `change` column from `calculate_vaf` function|
|`depth_vs_vaf_plot`|A scatterplot of depth vs VAF|A MAF formatted data.frame with the `change` column from `calculate_vaf` function|
|`snv_region_plot`|Genomic region ggplot barplot|An genomic region annotated MAF data.frame and has from `annotr_maf`|
|`cosmic_plot`|A violin plot of overlapping COSMIC vs non-COSMIC mutations|A MAF formatted data.frame with the `vaf` column from `calculate_vaf` function|
|`tmb_plot`|Plot Tumor Mutational Burden as a jitter plot|A data.frame with the tmb stats calculated from `calculate_tmb` function, `x_axis` to specify what variable to plot on the x-axis|
