# SNV caller comparison analysis

This analysis evaluates [MAF files](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from different SNV callers, compares their output, and creates a [consensus mutation file](./results/consensus/consensus_mutation.maf.tsv.zip).
This consensus mutation file is [MAF-like](#consensus-mutation-call) meaning it is TSV file that contains many of the fields of a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) but also some added calculations like [Variant Allele Fraction](#variant-allele-fraction-calculation) and some sample metadata information.

**Table of Contents**
* [How to run this pipeline](#how-to-run-this-pipeline)
* [General usage of scripts](#general-usage-of-scripts)
* [Individual caller evaluation](#individual-caller-evaluation)  
  * [Base change analysis](#base-change-analysis)  
  * [Variant allele fraction calculation](#variant-allele-fraction-calculation)  
  * [Genomic regional analyses](#genomic-regional-analyses)  
  * [Tumor mutation burden](#tumor-mutation-burden-calculation)
  * [COSMIC mutation overlap](#cosmic-mutation-overlap)  
* [Comparison analysis](#comparison-of-callers)
 * [Mutation IDs](#mutation-ids)
 * [Consensus mutation calls](#consensus-mutation-calls)
* [Overall file structure](#overall-file-structure)
* [Summary of functions](#summary-of-custom-functions)

## How to run the caller consensus analysis

To run the evaluations and comparisons of all the SNV callers, call the bash script:
```
bash run_caller_analysis.sh
```
This script will return results for each caller in the `plots` and `results` folder.
To see an overall summary report, look in the `results` folder for that caller.
(See [Overall File Structure](#overall-file-structure) for more details on
everything that is returned.)

The final results for the caller consensus are in the form of a notebook, `compare_snv_callers.nb.html`.
The consensus mutations themselves are saved to a [MAF-like file](#consensus-mutation-call) `consensus_mutation.maf.tsv.zip` to the `results/consensus`
results folder.
These are the mutations dubbed reliable enough to move forward with.

## General usage of scripts

**Overall notes about these scripts:**
- The scripts are sequential as noted by their number.
- All file path-related options assume the file path given is relative to `OpenPBTA-analysis`.
- By default, the scripts will not overwrite existing files of the same name. However,
this can be overridden with `--overwrite` option.

### 00-set_up.R

00-set_up.R creates the [annotation RDS file](#genomic-regional-analyses) and [COSMIC mutations file](#cosmic-mutation-overlap) that are used by
the subsequent scripts.
This set up script only needs to be run once and its three options are all relating to where the reference files should be stored.

**Option descriptions**
```
 --annot_rds : File path to where you would like the annotation_rds file to be
               stored
 --cosmic_og : Path to original COSMIC file. Can be .gz compressed. Will need to
               download this from COSMIC at https://cancer.sanger.ac.uk/cosmic/download
               These data are available if you register.
 --cosmic_clean : File path specifying where you would like the cleaned brain-related
                  COSMIC mutations file to be stored. This file is provided to you in
                  GitHub. Only coordinates are needed.
```

### 01-calculate_vaf_tmb.R

This script sets up the given MAF file and outputs three files ([VAF](#variant-allele-fraction-calculation),
[TMB](#tumor-mutation-burden-calculation), and [Regional analysis](#genomic-regional-analyses))
that are used to make an overall evaluation report in `02-run_eval.R`.

**Option descriptions**
```
 -label : Label to be used for folder and all output. eg. 'strelka2'. Default is 'maf'.
 -output : File path that specifies the folder where the output should go.
           New folder will be created if it doesn't exist.
 --file_format: What type of file format would you like the output as? Options are
               "rds" or "tsv". Default is "rds".
 --maf :  Relative file path to MAF file to be analyzed. Can be .gz compressed.
 --metadata : Relative file path to original metadata file.
 --annot_rds : Relative file path to annotation object RDS file to be analyzed.
 --bed_wgs : File path that specifies the caller-specific BED regions file.
 --bed_wxs : File path that specifies the WXS BED regions file.
 --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
 --no_region : If used, regional analysis will not be done.
```

### 02-run_eval.R

This script takes the files output by `01-calculate_vaf_tmb.R` and makes five
plots ([base_change](#base-change-analysis), [depth_vs_vaf](#variant-allele-fraction-calculation), [cosmic_plot](#cosmic-mutation-overlap), [snv_region](#genomic-regional-analyses), and [tmb_plot](#tumor-mutation-burden-calculation)) that are put into an overall report.

**Option descriptions**
```
 --label : Label to be used for folder and all output. eg. 'strelka2'. Optional.
           Default is 'maf'
 --plot_type : Specify what kind of plots you want printed out. Must be
               compatible with ggsave. eg pdf. Default is png
 --vaf : Folder from 01-calculate_vaf_tmb.R following files:
                                             <caller_name>_vaf.<file_format>
                                             <caller_name>_region.<file_format>
                                             <caller_name>_tmb.<file_format>
 --file_format: What type of file format were the vaf and tmb files saved as? Options are
                "rds" or "tsv". Default is "rds".
 --output : Where you would like the output from this script to be stored.
 --strategy : Specify whether you would like WXS and WGS separated for the plots.
              Analysis is still done on all data in the MAF file regardless.
              Acceptable options are 'wgs', 'wxs' or 'both', both for if you
              don't want to separate them. Default is both.
 --cosmic : Relative file path to COSMIC file to be analyzed.
 --overwrite : If TRUE, will overwrite any reports of the same name. Default is
              FALSE
  --no_region : If used, regional analysis will not be done.
```

# Individual Caller Evaluation

The first step in this analysis is an individual evaluation of each MAF from each caller.
This analysis prints out a report that is further subdivided by the WGS and WXS samples.

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
* `results/<caller_name>/<caller_name>_vaf.rds`
* `plots/<caller_name>/<caller_name>_<strategy>_base_change.png`

### Variant Allele Fraction Calculation

Calculate variant allele fraction (VAF) for each variant.
```
vaf = (t_alt_count) / (t_ref_count + t_alt_count)
```
This is following the [code used in
`maftools`](https://github.com/PoisonAlien/maftools/blob/1d0270e35c2e0f49309eba08b62343ac0db10560/R/plot_vaf.R#L39).
The VAF calculations and other special variables are added to the MAF fields and written to a TSV ending in `_vaf` in the caller's results folder.

 *Output for this analysis*
 * `results/<caller_name>/<caller_name>_vaf.rds`
 * `plots/<caller_name>/<caller_name>_<strategy>_depth_vs_vaf.png`

### Genomic Regional Analyses

To analyze what genomic regions the variants are from, I used [Annotatr
package](https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html) to obtain hg38 genome annotations.
This Annotatr object is stored as an RDS file: `hg38_genomic_region_annotations.rds` in the `scratch` directory.
Mutations are assigned all annotations that they overlap (using `GenomicRanges::overlap`).

*Output for this analysis*
* `results/<caller_name>/<caller_name>_regions.rds`
* `plots/<caller_name>/<caller_name>_<strategy>_snv_regions.png`

### Tumor Mutation Burden Calculation

To calculate TMB, the sum of the bases included in the WXS or WGS BED regions are used as the denominator, depending on the sample's processing strategy.

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
The sample-wise TMB calculations written to a TSV ending in `_tmb` in the caller's results folder.

*Output for this analysis*
* `results/<caller_name>/<caller_name>_tmb.rds`
* `plots/<caller_name>/<caller_name>_<strategy>_tmb_plot.png`

### COSMIC Mutation Overlap

The COSMIC mutation data were obtained from https://cancer.sanger.ac.uk/cosmic/download
*To run this analysis, you need to obtain these data.*
The full, unfiltered somatic mutations file `CosmicMutantExport.tsv` for grch38 is used here and the genomic coordinates is arranged to be in BED format.
The COSMIC set is filtered down to only mutations detected in brain-related
samples using the `Site subtype 1` field.
COSMIC mutations are overlapped with the present data's mutations using `GenomicRanges`.
The outcome of this overlap is added to the VAF data.frame with two `TRUE/FALSE` columns:
`overlap_w_cosmic` is TRUE for mutations that overlap with COSMIC mutations, while `same_as_cosmic` is TRUE when the base change summary is also identical.
The VAF for mutations that are or are not overlapping with COSMIC mutations are then plotted in a violin plot.

*Output for this analysis*
* `results/<caller_name>/<caller_name>_vaf.rds`
* `plots/<caller_name>/<caller_name>_<strategy>_cosmic_plot.png`

## Comparison of Callers

After running an initial evaluation and set up of each of the callers' MAF files,
the `compare_snv_callers.Rmd` notebook can be run to get final results on consensus
mutation calls.  

### Mutation IDs  

In order to compare mutations across callers, I created a `mutation_id` from combining information from standard MAF fields.
This was done in the `01-calculate_vaf_tmb.R` script using the `set_up_maf `
function.

`mutation_id` is a concatenation of:  
* `Hugo_Symbol`  
* [`change`](#base-change-analysis)  
* `Start_Position`  
* `Tumor_Sample_Barcode` (the sample ID)  

If mutation_id's are identical among MAF files, they are considered the same.

### Consensus mutation call

After the comparisons amongst the callers, VarDict proved to be too unreliable and called low VAF mutations.
Moving forward, mutations that were identified by Lancet, Mutect2 and Strelka2 were included in the final list of mutations for each sample.
The consensus mutations themselves are saved to a MAF-like file `consensus_mutation.maf.tsv.zip` to the `consensus`.
It is being called a "MAF-like" file because it has many of the same fields as a MAF file but..  
  - Does not contain the version string in the first row   
  - Has extraneous annotation data has been removed (columns with all `NA`s)  
  - Has VAF calculations and other variables that are calculated by the [`set_up_maf` function](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/util/wrangle_functions.R#L11).

## Overall file structure
```
OpenPBTA-analysis
├── analyses
│   └── snv-callers
│       ├── run_caller_analysis.sh
│       ├── compare_snv_callers.Rmd
│       ├── scripts
│       │   ├── 00-set-up.R
│       │   ├── 01-calculate_vaf_tmb.R
│       │   └── 02-run_eval.R
│       ├── util
│       │    ├── plot_functions.R
│       │    └── wrangle_functions.R
│       ├── results
│       │   ├── consensus
│       │   │   ├── consensus_mutation_tmb.tsv
│       │   │   └── consensus_mutation.maf.tsv
│       │   ├── lancet
│       │   │   ├── lancet_vaf.rds
│       │   │   ├── lancet_vaf.rds.zip
│       │   │   ├── lancet_tmb.rds
│       │   │   ├── lancet_tmb.rds.zip
│       │   │   ├── lancet_region.rds
│       │   │   ├── lancet_region.rds.zip
│       │   │   ├── lancet_wxs_report.html
│       │   │   ├── lancet_wxs_report.Rmd
│       │   │   ├── lancet_wgs_report.html
│       │   │   ├── lancet_wgs_report.Rmd
│       │   │   └── lancet_metadata_filtered.rds
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
│       ├── ref_files
│       │   ├── hg38_genomic_region_annotation.rds
│       │   └── brain_cosmic_variants_coordinates.tsv
│       └── template
│           ├── variant_caller_report_no_region_template.Rmd
│           └── variant_caller_report_template.Rmd
├── data
```

## Summary of custom functions

### Wrangling functions
|Function Name|Output created|Main Arguments|
|-------------|--------------|---------|
|`set_up_maf`|Columns: VAF, mutation_id, base_change, change| A MAF formatted data.frame|
|`maf_to_granges`|A `GenomicRanges` format object|A MAF formatted data.frame|
|`wxs_bed_filter`|A filtered MAF df with only mutations within the provided BED regions|A MAF formatted data.frame and a BED formatted data.frame|
|`calculate_tmb`|A data.frame with TMB stat per sample |A `wxs_bed_filter`ed MAF formatted data.frame, a WGS and WXS genome sizes|
|`annotr_maf`|A data.frame with genomic annotations for the provided MAF data.frame|A MAF formatted data.frame and a [built AnnotatR annotation object](https://rdrr.io/bioc/annotatr/man/build_annotations.html)|
|`find_cosmic_overlap`|Find overlap with [COSMIC](https://cancer.sanger.ac.uk/cosmic) mutations|A MAF formatted data.frame with `change` and variable from the `set_up_maf` function|

### Plotting functions

All plotting functions use `ggplot2` and all use a argument: `exp_strategy` that
determines whether to plot WGS or WXS samples only or both. Default is to plot
both.

|Function Name|Plot output|Main Arguments|
|-------------|-----------|---------|
|`base_change_plot`|Base change ggplot barplot|A MAF formatted data.frame with the `change` column from `set_up_maf` function|
|`depth_vs_vaf_plot`|A scatterplot of depth vs VAF|A MAF formatted data.frame with the `change` column from `set_up_maf` function|
|`snv_region_plot`|Genomic region ggplot barplot|An genomic region annotated MAF data.frame and has from `annotr_maf`|
|`cosmic_plot`|A violin plot of overlapping COSMIC vs non-COSMIC mutations|A MAF formatted data.frame with the `vaf` column from `set_up_maf` function|
|`tmb_plot`|Plot Tumor Mutational Burden as a jitter plot|A data.frame with the tmb stats calculated from `calculate_tmb` function, `x_axis` to specify what variable to plot on the x-axis|
