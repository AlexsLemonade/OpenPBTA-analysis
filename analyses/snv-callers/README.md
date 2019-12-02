# SNV caller comparison analysis

This analysis evaluates [MAF files](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from different SNV callers, compares their output, and creates a [consensus mutation file](./results/consensus/consensus_mutation.maf.tsv.zip).
This consensus mutation file is [MAF-like](#consensus-mutation-call) meaning it is TSV file that contains many of the fields of a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) but also some added calculations like [Variant Allele Fraction](#variant-allele-fraction-calculation) and some sample metadata information.

**Table of Contents**
* [How to run this pipeline](#how-to-run-this-pipeline)
* [General usage of scripts](#general-usage-of-scripts)
* [Individual caller evaluation](#individual-caller-evaluation)  
  * [Base change analysis](#base-change-analysis)  
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

### 01-setup_db.py

**Argument descriptions**
```
  -d DB_FILE, --db-file DB_FILE
     Path of the database file to use or create. Defaults to `data.sqlite`.
   --strelka-file STRELKA_FILE
     Path of the MAF formatted data file from the strelka2 caller(TSV).
   --mutect-file MUTECT_FILE
     Path of the MAF formatted data file from the mutect2 caller(TSV).
   --lancet-file LANCET_FILE
     Path of the MAF formatted data file from the lancet caller(TSV).
   --vardict-file VARDICT_FILE
     Path of the MAF formatted data file from the vardict caller(TSV).
   --meta-file META_FILE, --hist-file META_FILE
     Path of the metadata/histology data file(TSV).
   --overwrite           Overwrite tables that may already exist.
```

### 02-merge_callers.R

```
 --db_file : Path to sqlite database file made from 01-setup_db.py
 --output_file : File path and file name of where you would like the MAF-like
                 output from this script to be stored.
 --vaf_filter: Optional Variant Allele Fraction filter. Specify a number; any
               mutations with a VAF that are NA or below this number will be
               removed from the vaf data.frame before it is saved to a TSV file.
 --overwrite : If TRUE, will overwrite any reports of the same name. Default is
              FALSE
```
### 03-calculate_tmb.R

This script sets up the given MAF file and outputs three files ([VAF](#variant-allele-fraction-calculation),
[TMB](#tumor-mutation-burden-calculation), and [Regional analysis](#genomic-regional-analyses))
that are used to make an overall evaluation report in `02-run_eval.R`.

**Option descriptions**
```
 --consensus : File path to the MAF-like file.
 --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
 --bed_wgs : File path that specifies the caller-specific BED regions file.
             Assumes from top directory, 'OpenPBTA-analysis'.
 --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
             is given from top directory of 'OpenPBTA-analysis'
 --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
```
### Base change analysis

To evaluate base changes I summarized standard MAF fields as two new variables:
The `base_change` variable that indicates the exact change in bases from
concatenating `Reference_Allele`, `>`, and `Allele`.
The `change` variable is made from the `base_change` variable but groups
together deletions, insertions, and long (more than a SNV) as their own groups.

### Variant Allele Fraction Calculation

Calculate variant allele fraction (VAF) for each variant.
This is done in `03-calculate_tmb.R`.

```
vaf = (t_alt_count) / (t_ref_count + t_alt_count)
```
This is following the [code used in
`maftools`](https://github.com/PoisonAlien/maftools/blob/1d0270e35c2e0f49309eba08b62343ac0db10560/R/plot_vaf.R#L39).
The VAF calculations and other special variables are added to the MAF fields and written to a file ending in `_vaf` in the caller's results folder.

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
## Comparison of Callers

After running an initial evaluation and set up of each of the callers' MAF files,
the `compare_snv_callers.Rmd` notebook can be run to get final results on consensus
mutation calls.  

### Mutation IDs  

In order to compare mutations across callers, the data tables for each caller
were indexed by: `Chromosome` `Start_Position` `Reference_Allele` and `Allele`.
Meaning that across callers, if all these fields were identical, they were considered to be the same mutation.
This was done in the `01-setup_db.py` script using the function.

### Summary of consensus files:

- `consensus_mutation.maf.tsv` - Mutations that were called by all three of these callers for a given sample are saved to this file.
This file is [MAF-like](#consensus-mutation-call) meaning it is TSV file that contains many of the fields of a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) but also some added calculations like [Variant Allele Fraction](#variant-allele-fraction-calculation) and some sample metadata information.
These files combine the [MAF file data](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from 3 different SNV callers: [Mutect2](https://software.broadinstitute.org/cancer/cga/mutect), [Strelka2](https://github.com/Illumina/strelka), and [Lancet](https://github.com/nygenome/lancet).
See the methods on the callers' settings [here](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-single-nucleotide-variant-calling) and see the methods of this caller analysis and comparison [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers).

It is "MAF-like" file because it has many of the same fields as a MAF file but..  
  - Does not contain the version string in the first row   
  - Has extraneous annotation data has been removed (columns with all `NA`s)  
  - Has VAF calculations and other variables that are calculated by the [`set_up_maf` function](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/util/wrangle_functions.R#L11).

- `consensus_mutation_tmb.tsv` - After the consensus mutations were identified, Tumor mutation burden was recalculated for each sample from this mutation set.
See this [analysis' folder](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#tumor-mutation-burden-calculation) for details on these methods.

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
