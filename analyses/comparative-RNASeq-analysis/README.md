# Comparative RNA-Seq analysis

**Contents**

- [Purpose](#purpose)
- [Usage](#usage)
- [Limitations and Requirements](#limitations-and-requirements)
- [Future Updates](#future-updates)

## Purpose
The comparative-RNAseq-analysis module implements the outlier analysis workflow published in [Vaske et al. Jama Open Network. 2019](https://jamanetwork.com/journals/jamanetworkopen/article-abstract/2753519), which highlights genes within each sample whose expression is an outlier compared to the expression distribution of the dataset as a whole. This workflow:
  - Creates correlation matrices for polyA samples and ribodeplete (stranded) samples using gene expression data.
  - Generates gene outlier threshold values for ribodeplete samples.
  - Lists outlier genes for each ribodeplete sample.

## Usage
The input file must be an RNA-Seq TPM gene expression matrix in the [RDS](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) file format. Currently available matrices that meet this format are `pbta-gene-expression-rsem-tpm.polya.rds` and `pbta-gene-expression-rsem-tpm.stranded.rds`.

Command lines follow, using as example the stranded (ribodeplete) dataset.
Available flags:
  - `--verbose` toggles verbose output
  - `--output-prefix MyDataset` prepends MyDataset to output filenames to identify runs. Subsequent steps in the same run must use the same prefix.
  - `--scratch ../../scratch` provides path to scratch dir where intermediate files shared between steps can be read and written.
  - `--results ./results` provides path to final results dir.

### 01 - Correlation Matrix
Generates the correlation matrix and filtered gene list.

```
./scripts/run_in_ci.sh \
  python3 analyses/comparative-RNASeq-analysis/01-correlation-matrix.py \
    ../../data/pbta-gene-expression-rsem-tpm.stranded.rds \
    --output-prefix rsem-tpm-stranded- \
    --verbose
```

Input file:
```
data/pbta-gene-expression-rsem-tpm.stranded.rds
```

Output files:
```
scratch/rsem-tpm-stranded-all_by_all_correlations.rds
scratch/rsem-tpm-stranded-filtered_genes_to_keep.rds
scratch/rsem-tpm-stranded-log2-normalized.rds
```

### 02 -  Thresholds and Outliers
Generates outlier thresholds and matrix of outlier genes.

```
./scripts/run_in_ci.sh \
  python3 analyses/comparative-RNASeq-analysis/02-thresholds-and-outliers.py \
    --prefix rsem-tpm-stranded- \
    --verbose
```

Additional flags for this step:
  - `--iqr-multiplier 1.5` sets the interquartile range multiplier for the Tukey outlier calculation


Input files (detected via command-line prefix provided):
```
scratch/rsem-tpm-stranded-log2-normalized.rds
scratch/rsem-tpm-stranded-filtered_genes_to_keep.rds
```

Output files:
```
results/rsem-tpm-stranded-gene_expression_outliers.tsv.gz
```

## Limitations and requirements
Because the per-sample results of this analysis are dependent on the entire dataset, all samples in the dataset must meet certain standards for the outliers to be meaningful. *(Currently, these standards are not being enforced.)*
  - All samples must pass a quality control check.
  - Dataset must contain only tumor samples; no normal, cell line, etc data.
  - All samples in the dataset must have the same library preparation for their gene expression to be comparable. (Eg, polyA selection, ribodepletion, or hybrid capture).

### Software dependencies
The analysis uses python 3 and requires the following libraries. Version numbers
are those currently in use and earlier or later versions may also be acceptable but have not been tested.
```
numpy (1.17.3)
pandas (0.25.3)
scipy (1.3.2)
scikit-learn (0.19.1)
pyreadr (0.2.1)
```
In addition, a `utils` module is included for certain shared functions.

## Future updates
- The next iteration of this module will use the `pbta-histologies.tsv` file to filter the input datasets to only tumor samples. It will also filter samples to those which meet a to-be-described quality control standard.

