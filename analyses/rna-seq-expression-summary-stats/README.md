## Calculate TPM summary statistics for each cancer group and cohort

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

For each cancer group and cohort, calculate TPM means, standard deviations, z-scores, and ranks.

### Methods

Select independent RNA-seq samples using `independent-specimens.rnaseq.primary.eachcohort.tsv` in the results of the `independent-samples` analysis module.

Group all samples with one of the following two methods:

- Each `sample_group` includes all samples in one `cancer_group` and one `cohort`. Call this grouping method `each-cohort`.
- Each `sample_group` includes all samples in one `cancer_group` and all `cohort`(s). Call this grouping method `all-cohorts` (or `PedOT`).

For each grouping method, compute summary statistics as following:

For all `sample_group`s with >= 5 samples, generate the following summary statistics:

- TPM means of each gene across all selected samples, and denote this vector as `mean_TPM_vector`.
- TPM standard deviations of each gene across all selected samples.
- z-scores of each gene across all genes, computed as `z_score_vector = (mean_TPM_vector - mean(mean_TPM_vector)) / sd(mean_TPM_vector)`. Call these z-scores as `mean_tpm_cancer_group_wise_z_scores` in the filenames.
- Ranks of mean TPM of each gene across all selected samples, which takes one of the following four values: `Highest expressed 25%`, `Expression between upper quartile and median`, `Expression between median and lower quartile`, and `Lowest expressed 25%`. If multiple genes have the same mean TPM value, their tied rank is the lowest rank, in order to be conservative on the description of their expression levels.

Combine each type of the summary statistics vectors into a table, with rows as genes, and columns as `sample_group`s.

Denote the (`n_genes`, `n_sample_groups`) `mean_TPM_vector` combined matrix as `mean_TPM_matrix`.

Generate z-scores across all `sample_groups` as `z_score_matrix = (mean_TPM_matrix - rowMeans(mean_TPM_matrix)) / rowSD(mean_TPM_matrix)`. Call these z-scores as `mean_tpm_gene_wise_z_scores` in the filenames.

Call aforementioned (`n_genes`, `n_sample_groups`) tables as wide tables.

Generate long summary statistic tables for converting to JSON format, with each row as a tab-delimited record of the following columns, as suggested by @jharenza and @taylordm at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/27#issuecomment-868035006>.

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

### Results

The `NA`/`NaN`s in result tables are represented with blank string `''`s.

#### `all-cohorts`(/`PedOT`) summary statistics tables

The following wide tables are generated using the methods described above. Rows are genes. Columns are `sample_group`s, except that the first two columns are gene symbol and gene Ensembl ID.

In `all-cohorts` result tables, a `sample_group` is a string that concatenates a `cancer_group` and `___all_cohorts`, e.g. `Meningioma___all_cohorts`, `Neuroblastoma___all_cohorts`, and `Diffuse midline glioma___all_cohorts`.


If one gene symbol matches to multiple Ensembl IDs, the value of the Ensembl ID column is a comma separated list of all Ensembl IDs, e.g. `ENSG00000206952,ENSG00000281910`. In `inpiut/ens_symbol.tsv`, `SNORA50A` is mapped to both `ENSG00000206952` and `ENSG00000281910`.

- `results/cancer_group_all_cohort_mean_tpm.tsv`
- `results/cancer_group_all_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv`

The samples used in each `sample_group` are listed in `results/cancer_group_all_cohort_sample_metadata.tsv`. The columns are as following:

- `sample_group`
- included one or more `cohort`s in the `sample_group`
- the number of samples in the `sample_group`
- the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `sample_group`

#### `each-cohort` summary statistics tables

The following wide tables are generated using the methods described above. Rows are genes. Columns are `sample_group`s, except that the first two columns are gene symbol and gene Ensembl ID.

In `each-cohort` result tables, a `sample_group` is a string that concatenates a `cancer_group` and a `cohort` by `___`, e.g. `Meningioma___PBTA`, `Neuroblastoma___GMKF`, and `Diffuse midline glioma___PBTA`.

If one gene symbol matches to multiple Ensembl IDs, the value of the Ensembl ID column is a comma separated list of all Ensembl IDs, e.g. `ENSG00000206952,ENSG00000281910`. In `inpiut/ens_symbol.tsv`, `SNORA50A` is mapped to both `ENSG00000206952` and `ENSG00000281910`.

- `results/cancer_group_individual_cohort_mean_tpm.tsv`
- `results/cancer_group_individual_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_gene_wise_z_scores.tsv`

The samples used in each `sample_group` are listed in `results/cancer_group_individual_cohort_sample_metadata.tsv`. The columns are as following:

- `sample_group`
- the number of samples in the `sample_group`
- the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `sample_group`

#### Long summary statistics tables

The following long tables are generated using the wide tables in TSV and JSON formats.

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

If the record is generated using the `all-cohorts`(/`PedOT`) grouping method, the `cohort` column takes the value of `all_cohorts`.

Different from the wide tables, the `gene_id` does not contain comma-separated ENSG IDs. If one gene symbol matches to multiple Ensembl IDs, each Ensembl gene ID will become one row in the long table. For example:

- Wide table:

```text
# cancer_group_individual_cohort_mean_tpm.tsv
gene_symbol    gene_id    Adamantimomatous craniopharyngioma___PBTA    Atypical meningioma___PBTA    Atypical Teratoid Rhabdoid Tumor___PBTA
CDR1    ENSG00000184258,ENSG00000281508    111.1765    31.932000000000002    332.2896153846154
```

- Long table:

```text
# long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz
CDR1    ENSG00000184258    Adamantimomatous craniopharyngioma    PBTA    111.1765    181.70448082054318    0.11408800038186823    Highest expressed 25%    20
CDR1    ENSG00000281508    Adamantimomatous craniopharyngioma    PBTA    111.1765    181.70448082054318    0.11408800038186823    Highest expressed 25%    20
CDR1    ENSG00000184258    Atypical meningioma    PBTA    31.932000000000002    50.10366423725914    0.016467672549703015    Highest expressed 25%    5
CDR1    ENSG00000281508    Atypical meningioma    PBTA    31.932000000000002    50.10366423725914    0.016467672549703015    Highest expressed 25%    5
CDR1    ENSG00000184258    Atypical Teratoid Rhabdoid Tumor    PBTA    332.2896153846154    1075.378786603049    0.37106139758788453    Highest expressed 25%    26
CDR1    ENSG00000281508    Atypical Teratoid Rhabdoid Tumor    PBTA    332.2896153846154    1075.378786603049    0.37106139758788453    Highest expressed 25%    26
```

### Usage

1. Change working directory to local `OpenPBTA-analysis`.
2. Download data using `bash download-data.sh`.
3. Run this analysis module in the continuous integration (CI) docker image using `./scripts/run_in_ci.sh bash analyses/rna-seq-expression-summary-stats/run-rna-seq-expression-summary-stats.sh`.

### Module structure

```text
.
├── 01-tpm-summary-stats.R
├── README.md
├── results
│   ├── cancer_group_all_cohort_mean_tpm.tsv
│   ├── cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv
│   ├── cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv
│   ├── cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv
│   ├── cancer_group_all_cohort_sample_metadata.tsv
│   ├── cancer_group_all_cohort_standard_deviation_tpm.tsv
│   ├── cancer_group_individual_cohort_mean_tpm.tsv
│   ├── cancer_group_individual_cohort_mean_tpm_cancer_group_wise_quantiles.tsv
│   ├── cancer_group_individual_cohort_mean_tpm_cancer_group_wise_z_scores.tsv
│   ├── cancer_group_individual_cohort_mean_tpm_gene_wise_z_scores.tsv
│   ├── cancer_group_individual_cohort_sample_metadata.tsv
│   ├── cancer_group_individual_cohort_standard_deviation_tpm.tsv
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.json.gz
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.jsonl.gz
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz
│   ├── long_n_tpm_mean_sd_quantile_group_wise_zscore.json.gz
│   ├── long_n_tpm_mean_sd_quantile_group_wise_zscore.jsonl.gz
│   └── long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz
└── run-rna-seq-expression-summary-stats.sh
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
- `../../data/ensg-hugo-rmtl-v1-mapping.tsv`
- `../../data/efo-mondo-map.tsv`
- `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`

Output:

- `results/cancer_group_all_cohort_sample_metadata.tsv`
- `results/cancer_group_all_cohort_mean_tpm.tsv`
- `results/cancer_group_all_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv`

- `results/cancer_group_individual_cohort_sample_metadata.tsv`
- `results/cancer_group_individual_cohort_mean_tpm.tsv`
- `results/cancer_group_individual_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_gene_wise_z_scores.tsv`

- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.json`
- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv`
- `long_n_tpm_mean_sd_quantile_group_wise_zscore.json`
- `long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv`
