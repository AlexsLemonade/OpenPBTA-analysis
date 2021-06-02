# Comparative RNA-Seq analysis

**Contents**

- [Purpose](#purpose)
- [Usage](#usage)
  - [01 - Correlation Matrix](#01---correlation-matrix)
  - [02 -  Thresholds and Outliers](#02----thresholds-and-outliers)
- [Limitations and requirements](#limitations-and-requirements)
  - [Software dependencies](#software-dependencies)
- [Future updates](#future-updates)

## Purpose
The comparative-RNAseq-analysis module implements the outlier analysis workflow published in [Vaske et al. Jama Open Network. 2019](https://jamanetwork.com/journals/jamanetworkopen/article-abstract/2753519), which highlights genes within each sample whose expression is an outlier compared to the expression distribution of the dataset as a whole. This workflow:
  - Creates correlation matrices for polyA samples and ribodeplete (stranded) samples using gene expression data.
  - Generates gene outlier threshold values for ribodeplete samples.
  - Lists outlier genes for each ribodeplete sample.

## Usage
The input file must be an RNA-Seq TPM gene expression matrix in the [RDS](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) file format. Currently available matrices that meet this format are `pbta-gene-expression-rsem-tpm.polya.rds` and `pbta-gene-expression-rsem-tpm.stranded.rds`.

Command lines follow, using as example the stranded (ribodeplete) dataset.

### 01 - Correlation Matrix
Filters input samples to only tumor samples which are MEND QC pass using the histology spreadsheet
and qc manifest and results.
Generates the correlation matrix and filtered gene list.

```
python3 01-correlation-matrix.py \
  ../../data/pbta-gene-expression-rsem-tpm.stranded.rds \
  --clinical-path ../../data/pbta-histologies.tsv \
  --qc-manifest-path ../../data/pbta-mend-qc-manifest.tsv \
  --qc-results-path ../../data/pbta-mend-qc-results.tar.gz \
  --prefix rsem-tpm-stranded- \
  --verbose
```

Required flags:
  - First argument must be the path to the input expression .Rds file
  - `--scratch ../../scratch` provides path to scratch dir where intermediate files shared between steps can be read and written.

Optional flags:
  - `--verbose` enables verbose output.
  - `--prefix MyDataset` prepends MyDataset to output filenames to identify runs. Subsequent steps in the same run must use the same prefix.
  - `--nofilter` causes the sample filtering step to be skipped.
If the sample filtering step is *not* skipped, these flags are required:
  - `--clinical-path` path to clinical tsv file with tumor status columns, eg `pbta-histologies.tsv`.
  - `--qc-manifest-path` path to the qc manifest mapping sample IDs to filenames, eg `pbta-mend-qc-manifest.tsv`.
  - `--qc-results-path` path to the tarball of  qc results files, eg `pbta-mend-qc-results.tar.gz`.


Input files:
```
data/pbta-gene-expression-rsem-tpm.stranded.rds
data/pbta-histologies.tsv
data/pbta-mend-qc-manifest.tsv
data/pbta-mend-qc-results.tar.gz
```

Output files:
```
scratch/rsem-tpm-stranded-all_by_all_correlations.feather
scratch/rsem-tpm-stranded-gene-filters.feather
scratch/rsem-tpm-stranded-filtered-log2-normalized.feather
```

### 02 -  Thresholds and Outliers
Generates outlier thresholds and matrix of outlier genes.

```
python3 02-thresholds-and-outliers.py \
  --prefix rsem-tpm-stranded- \
  --results results \
  --verbose
```

Required flags:
 - `--scratch ../../scratch` provides path to scratch dir where intermediate files shared between steps can be read and written.
 - `--results results` provides path to final results dir.

Optional flags:
  - `--verbose` enables verbose output
  - `--prefix MyDataset` prepends MyDataset to input and output filenames to identify runs. Subsequent steps in the same run must use the same prefix.
  - `--iqr-multiplier 1.5` sets the interquartile range multiplier for the Tukey outlier calculation

Input files (detected via the provided prefix):
```
scratch/rsem-tpm-stranded-filtered-log2-normalized.feather
scratch/rsem-tpm-stranded-gene-filters.feather
```

Output files:
```
results/rsem-tpm-stranded-gene_expression_outliers.tsv.gz
```

## Limitations and requirements
Because the per-sample results of this analysis are dependent on the entire dataset, all samples in the dataset must meet certain standards for the outliers to be meaningful:
  - All samples must pass the MEND quality control check, which confirms that the sample has at least 10 million Mapped Exonic Non-Duplicate reads.
  - Dataset must contain only tumor samples; no normal, cell line, etc data. The filter applied via the `pbta-histologies.tsv` file to enforce this is:
```
sample_type == Tumor
composition == Solid Tissue
experimental_strategy == RNA-Seq
```
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
pyarrow (0.16.0)
```
In addition, a `utils` module is included for certain shared functions.

## Future updates
- Trends in outlier results
- PolyA vs RiboDeplete classifier
