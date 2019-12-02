# SNV caller comparison analysis

This analysis evaluates [MAF files](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from different SNV callers, compares their output, and creates a [consensus mutation file](./results/consensus/consensus_mutation.maf.tsv.zip).
This consensus mutation file is [MAF-like](#consensus-mutation-call) meaning it is TSV file that contains many of the fields of a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) but also has [VAF](#variant-allele-fraction-calculation)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [How to run the caller consensus analysis](#how-to-run-the-caller-consensus-analysis)
  - [Summary of consensus files:](#summary-of-consensus-files)
- [Summary of Methods](#summary-of-methods)
  - [Mutation Comparisons](#mutation-comparisons)
  - [Variant Allele Fraction Calculation](#variant-allele-fraction-calculation)
  - [Tumor Mutation Burden Calculation](#tumor-mutation-burden-calculation)
- [General usage of scripts](#general-usage-of-scripts)
  - [01-setup_db.py](#01-setup_dbpy)
  - [02-merge_callers.R](#02-merge_callersr)
  - [03-calculate_tmb.R](#03-calculate_tmbr)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

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

## Summary of Methods

### Mutation Comparisons

In order to compare mutations across callers, the data tables for each caller
were indexed by: `Chromosome` `Start_Position` `Reference_Allele` and `Allele`.
Meaning that across callers, if all these fields were identical, they were considered to be the same mutation.
This was done in the `01-setup_db.py` script using the function.

### Variant Allele Fraction Calculation

Calculate variant allele fraction (VAF) for each variant.
This is done in `01-setup_db.py`.

```
vaf = (t_alt_count) / (t_ref_count + t_alt_count)
```
This is following the [code used in
`maftools`](https://github.com/PoisonAlien/maftools/blob/1d0270e35c2e0f49309eba08b62343ac0db10560/R/plot_vaf.R#L39).

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

## General usage of scripts

**Overall notes about these scripts:**
- The scripts are sequential as noted by their number.
- All file path-related options assume the file path given is relative to `OpenPBTA-analysis`.
- By default, the scripts will not overwrite existing files of the same name. However,
this can be overridden with `--overwrite` option.

### 01-setup_db.py

Creates and/or fills a database of variant calls that will be used by subsequent calls.
Note: requires `pandas` to be installed, and expects python3
All arguments are optional; only the included tables will be affected.

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

Using the database created by `01-setup_db.py`, merge callers' data files into consensus MAF-like file.

**Argument descriptions**
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

Using the consensus file created in `02-merge_callers.R`, calculate TMB for all
WGS and WXS samples.

**Argument descriptions**
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
