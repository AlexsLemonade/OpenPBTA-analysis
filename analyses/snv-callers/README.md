# SNV caller comparison analysis

This analysis evaluates [MAF files](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from different SNV callers, compares their output, and creates a [consensus mutation file](./results/consensus/consensus_mutation.maf.tsv.zip).
This consensus mutation file is [MAF-like](#consensus-mutation-call) meaning it is TSV file that contains the fields of a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) but adds [VAF](#variant-allele-fraction-calculation), but does not contain a starting comment line with a version number.

See the comparison results plots [here](https://cansavvy.github.io/openpbta-notebook-concept/snv-callers/compare_snv_callers_plots.nb.html).

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
This bash script will return:

- Comparison plots in a notebook: [`compare_snv_callers_plots.nb.html`](https://cansavvy.github.io/openpbta-notebook-concept/snv-callers/compare_snv_callers_plots.nb.html).
- A zip file containing:
  - `consensus_snv.maf.tsv` - is  [MAF-like file](#consensus-mutation-call) that contains the snvs that were called by all three of these callers for a given sample are saved to this file.
  These files combine the [MAF file data](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from 3 different SNV callers: [Mutect2](https://software.broadinstitute.org/cancer/cga/mutect), [Strelka2](https://github.com/Illumina/strelka), and [Lancet](https://github.com/nygenome/lancet).
  See the methods on the callers' settings [here](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-single-nucleotide-variant-calling) and see [the methods of this caller analysis and comparison below](#summary-of-methods).  
  - `consensus_snv_tmb_coding_only.tsv` - Tumor Mutation burden calculations using *coding only* mutations use the consensus of Lancet, Mutect2, and Strelka2. 
  - `consensus_snv_tmb_all.tsv` - Tumor Mutation burden calculations using *all* mutations use the consensus of Mutect2, and Strelka2. (Lancet was excluded because it has a [coding region bias in the way it was ran](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling)).

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

Creates and/or fills an SQLite database of variant calls that will be used by subsequent steps to find consensus mutations.
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

Using the database created by `01-setup_db.py`, merge callers' data files into consensus [MAF-like file](#snv-caller-comparison-analysis).

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
Two TMB files are created, one including *all snv* called by Strelka2 and Mutect2 (Lancet is excluded from this TMB calculation consensus because of a [coding region bias in the way it was ran](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling)), and a *coding snvs only* TMB calculation. 

**Argument descriptions**
```
 --consensus : File path to the MAF-like file.
 --db_file : Path to sqlite database file made from 01-setup_db.py
 --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
 --bed_wgs : File path that specifies the caller-specific BED regions file.
             Assumes from top directory, 'OpenPBTA-analysis'.
 --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
             is given from top directory of 'OpenPBTA-analysis'
 --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
```
