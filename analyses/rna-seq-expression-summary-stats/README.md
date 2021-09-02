## Calculate TPM summary statistics for each cancer group and cohort

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

- [Calculate TPM summary statistics for each cancer group and cohort](#calculate-tpm-summary-statistics-for-each-cancer-group-and-cohort)
  - [Purpose](#purpose)
  - [Methods](#methods)
  - [Results](#results)
    - [`all-cohorts` sample metadata table](#all-cohorts-sample-metadata-table)
    - [`each-cohort` sample metadata table](#each-cohort-sample-metadata-table)
    - [Long summary statistics tables](#long-summary-statistics-tables)
  - [Usage](#usage)
  - [Module structure](#module-structure)
  - [Analysis scripts](#analysis-scripts)
    - [01-tpm-summary-stats.R](#01-tpm-summary-statsr)
      - [Unit testing](#unit-testing)

### Purpose

For each cancer group and cohort, calculate TPM means, standard deviations, z-scores, and ranks.

### Methods

Select independent RNA-seq samples using `independent-specimens.rnaseq.primary.eachcohort.tsv` for `all cohorts` and `independent-specimens.rnaseq.primary.eachcohort.tsv` for `each cohort` in the results of the `independent-samples` analysis module.

Group all samples with one of the following two methods:

- Each `sample_group` includes all samples in one `cancer_group` and one `cohort`. Call this grouping method `each-cohort`.
- Each `sample_group` includes all samples in one `cancer_group` and all `cohort`(s). Call this grouping method `all-cohorts`.

For each grouping method, compute summary statistics as following:

For all `sample_group`s with >= 3 samples, generate the following summary statistics:

- TPM means of each gene across all selected samples, and denote this vector as `mean_TPM_vector`.
- TPM standard deviations of each gene across all selected samples.
- z-scores of each gene across all genes, computed as `z_score_vector = (mean_TPM_vector - mean(mean_TPM_vector)) / sd(mean_TPM_vector)`. Call these z-scores as `mean_tpm_cancer_group_wise_z_scores` in the filenames.
- Ranks of mean TPM of each gene across all selected samples, which takes one of the following four values: `Highest expressed 25%`, `Expression between upper quartile and median`, `Expression between median and lower quartile`, and `Lowest expressed 25%`. If multiple genes have the same mean TPM value, their tied rank is the lowest rank, in order to be conservative on the description of their expression levels.

Combine each type of the summary statistics vectors into a table, with rows as genes, and columns as `sample_group`s.

Denote the (`n_genes`, `n_sample_groups`) `mean_TPM_vector` combined matrix as `mean_TPM_matrix`.

Generate z-scores across all `sample_groups` as `z_score_matrix = (mean_TPM_matrix - rowMeans(mean_TPM_matrix)) / rowSD(mean_TPM_matrix)`. Call these z-scores as `mean_tpm_gene_wise_z_scores` in the filenames.

Generate long summary statistic tables for converting to [JSON Lines (JSONL)](<https://jsonlines.org/>) format, with each row as a tab-delimited record of the following columns, as suggested by @jharenza and @taylordm at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/27#issuecomment-868035006>.

- `gene_symbol`
- `RMTL`
- `gene_id`
- `cancer_group`
- `EFO`
- `MONDO`
- `n_samples`
- `cohort`
- `tpm_mean`
- `tpm_sd`
- `tpm_mean_cancer_group_wise_zscore`/`tpm_mean_gene_wise_zscore`
- `tpm_mean_cancer_group_wise_quantiles`

Output the long summary statistic tables in TSV and JSON formats using `readr::write_tsv()` and `jsonlite::write_json()` respectively.

Convert JSON files to JSONL format with [`jq`](<https://stedolan.github.io/jq/>).

### Results

The results are generated using v7 data release and independent sample lists.

The `NA`/`NaN`s in result tables are replaced with blank string `''`s.

#### `all-cohorts` sample metadata table

The samples used in each `all-cohorts` `sample_group` are listed in `results/cancer_group_all_cohorts_sample_metadata.tsv`. The columns are as following:

- `sample_group`
- included one or more `cohort`s in the `sample_group`
- the number of samples in the `sample_group`
- the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `sample_group`

#### `each-cohort` sample metadata table

The samples used in each `each-cohort` `sample_group` are listed in `results/cancer_group_each_cohort_sample_metadata.tsv`. The columns are as following:

- `sample_group`
- the number of samples in the `sample_group`
- the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `sample_group`

#### Long summary statistics tables

The following long tables are generated using the computed summary statistics in TSV and JSONL formats.

- `results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.json.gz`
- `results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz`
- `results/long_n_tpm_mean_sd_quantile_group_wise_zscore.json.gz`
- `results/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz`

Each row of the long table is a tab-delimited record of the following columns, as suggested by @jharenza and @taylordm at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/27#issuecomment-868035006>.

- `gene_symbol`
- `RMTL`
- `gene_id`
- `cancer_group`
- `EFO`
- `MONDO`
- `n_samples`
- `cohort`
- `tpm_mean`
- `tpm_sd`
- `tpm_mean_cancer_group_wise_zscore`/`tpm_mean_gene_wise_zscore`
- `tpm_mean_cancer_group_wise_quantiles`

If the record is generated using the `all-cohorts` grouping method, the `cohort` column takes the value of `all_cohorts`.

If one gene symbol matches to multiple Ensembl gene IDs, each mapped Ensembl gene ID will become one row in the long table. For example:

- Gene symbol CDR1 is mapped to the following two Ensemble gene IDs, ENSG00000184258 and ENSG00000281508.

- Long table rows:

```text
$ gunzip -c results/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz | grep -P '^CDR1\t.*(Atypical Teratoid Rhabdoid Tumor|Meningioma)' 
CDR1            ENSG00000184258 Atypical Teratoid Rhabdoid Tumor        EFO_1002008     MONDO_0020560   24      all_cohorts     356.4175        1117.43530723843        0.402923612719455  Highest expressed 25%
CDR1            ENSG00000281508 Atypical Teratoid Rhabdoid Tumor        EFO_1002008     MONDO_0020560   24      all_cohorts     356.4175        1117.43530723843        0.402923612719455  Highest expressed 25%
CDR1            ENSG00000184258 Meningioma      EFO_0003851     MONDO_0016642   14      all_cohorts     63.1678571428571        177.200539857598        0.0494325433009981      Highest expressed 25%
CDR1            ENSG00000281508 Meningioma      EFO_0003851     MONDO_0016642   14      all_cohorts     63.1678571428571        177.200539857598        0.0494325433009981      Highest expressed 25%
CDR1            ENSG00000184258 Atypical Teratoid Rhabdoid Tumor        EFO_1002008     MONDO_0020560   24      PBTA    356.4175        1117.43530723843        0.402923612719455       Highest expressed 25%
CDR1            ENSG00000281508 Atypical Teratoid Rhabdoid Tumor        EFO_1002008     MONDO_0020560   24      PBTA    356.4175        1117.43530723843        0.402923612719455       Highest expressed 25%
CDR1            ENSG00000184258 Meningioma      EFO_0003851     MONDO_0016642   14      PBTA    63.1678571428571        177.200539857598        0.0494325433009981      Highest expressed 25%
CDR1            ENSG00000281508 Meningioma      EFO_0003851     MONDO_0016642   14      PBTA    63.1678571428571        177.200539857598        0.0494325433009981      Highest expressed 25%
```

### Usage

1. Change working directory to local `OpenPedCan-analysis`.
2. Download data using `bash download-data.sh`.
3. Run this analysis module in the continuous integration (CI) docker image using `./scripts/run_in_ci.sh bash analyses/rna-seq-expression-summary-stats/run-rna-seq-expression-summary-stats.sh`.

### Module structure

```text
.
├── 01-tpm-summary-stats.R
├── README.md
├── results
│   ├── cancer_group_all_cohorts_sample_metadata.tsv
│   ├── cancer_group_each_cohort_sample_metadata.tsv
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.jsonl.gz
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz
│   ├── long_n_tpm_mean_sd_quantile_group_wise_zscore.jsonl.gz
│   └── long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz
├── run-rna-seq-expression-summary-stats.sh
├── run-tests.sh
└── tests
    ├── helper_import_function.R
    ├── test_data
    │   ├── test_import_function_empty.R
    │   └── test_import_function_non_empty.R
    ├── test_format_cohort_sample_counts.R
    └── test_helper_import_function.R
```

### Analysis scripts

#### 01-tpm-summary-stats.R

Usage:

```bash
Rscript --vanilla '01-tpm-summary-stats.R'
```

Input:

- `../../data/gene-expression-rsem-tpm-collapsed.rds`
- `../../data/histologies.tsv`
- `../../data/ensg-hugo-rmtl-mapping.tsv`
- `../../data/independent-specimens.rnaseq.primary.tsv`
- `../../data/independent-specimens.rnaseq.primary.eachcohort.tsv`

Output:

- `results/cancer_group_all_cohorts_sample_metadata.tsv`
- `results/cancer_group_each_cohort_sample_metadata.tsv`
- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.json`
- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv`
- `long_n_tpm_mean_sd_quantile_group_wise_zscore.json`
- `long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv`

##### Unit testing

The unit testing is implemented using the [`testthat`](https://testthat.r-lib.org/index.html) package version 2.1.1, as suggested by @jharenza and @NHJohnson in the reviews of PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/55>.

To run all unit tests, run `bash run-tests.sh` in the Docker image/container from any working directory. Following is an example run.

```text
$ bash analyses/rna-seq-expression-summary-stats/run-tests.sh
✔ |  OK F W S | Context
✔ |  15       | tests/test_format_cohort_sample_counts.R [0.3 s]
✔ |  21       | tests/test_helper_import_function.R

══ Results ═════════════════
Duration: 0.4 s

OK:       36
Failed:   0
Warnings: 0
Skipped:  0
Done running run-tests.sh
```

To add more tests, create additional `test*R` files under the `tests` directory, with available `test*R` files as reference.

Notes on the `testthat` unit testing framework:

- `testthat::test_dir("tests")` finds all `test*R` files under the `tests` directory to run, which is used in `run-tests.sh`.
- `testthat::test_dir("tests")` also finds and runs all `helper*R` files under the `tests` directory before running the `test*R` files.
- The working directory is `tests` when running the `helper*R` and `test*R` files through `testthat::test_dir("tests")`.
- In order to import a funciton for testing from an R file without running the whole file, a helper function `import_function` is defined at `tests/helper_import_function.R`, and the `import_function` is also tested in the `tests/test_helper_import_function.R` file.
- Even though the `testthat` 2.1.1 documentation of the `filter` parameter of `test_dir` function says that "Matching is performed on the file name after it's stripped of "test-" and ".R", the R code uses the following. Therefore, naming test files with `test_some_test_file.R` can be found by the `test_dir` function.
  - `"^test.*\\.[rR]$"` for finding test files in `find_test_scripts`
  - `sub("^test-?", "", test_names)`, `sub("\\.[rR]$", "", test_names)`, and `grepl(filter, test_names, ...)` for filtering test files in `testthat:::filter_test_scripts`.
