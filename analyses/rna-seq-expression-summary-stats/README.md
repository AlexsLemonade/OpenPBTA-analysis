## Calculate TPM summary statistics within each cancer group and cohort

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

Within each cancer group and cohort, calculate TPM means, standard deviations, z-scores, and ranks.

### Methods

Select independent RNA-seq samples using `independent-specimens.rnaseq.primary.tsv` in the results of the `independent-samples` analysis module.


For each `cancer_group`, select one of the following two sets of samples:

- Samples from all `cohort`s, e.g. PBTA and GMKF.
- Samples from each individual `cohort`.

If >= 5 samples are selected, generate the following summary statistics:

- TPM means of each gene across all selected samples, and denote this vector as `mean_TPM_vector`.
- TPM standard deviations of each gene across all selected samples.
- z-scores of each gene across all genes, computed as `z_score_vector = (mean_TPM_vector - mean(mean_TPM_vector)) / sd(mean_TPM_vector)`. Call these z-scores as `mean_tpm_cancer_group_wise_z_scores` in the filenames.
- Ranks of mean TPM of each gene across all selected samples, which takes one of the following four values: `Highest expressed 25%`, `Expression between upper quartile and median`, `Expression between median and lower quartile`, and `Lowest expressed 25%`. If multiple genes have the same mean TPM value, their tied rank is the lowest rank, in order to be conservative on the description of their expression levels.

Combine each type of the summary statistics vectors into a table, with rows as genes, and columns as `cancer_group_cohort`.

Denote the (`n_genes`, `n_cancer_groups`/`n_cancer_group_cohorts`) `mean_TPM_vector` combined matrix as `mean_TPM_matrix`.

Generate z-scores across all `cancer_groups`/`cancer_group_cohorts` as `z_score_matrix = (mean_TPM_matrix - rowMeans(mean_TPM_matrix)) / rowSD(mean_TPM_matrix)`. Call these z-scores as `mean_tpm_gene_wise_z_scores` in the filenames.

Call aforementioned (`n_genes`, `n_cancer_groups`/`n_cancer_group_cohorts`) tables as wide tables.

Generate long summary statistic tables for converting to JSON format, with each row as a tab-delimited record of the following columns, as suggested by @jharenza and @taylordm at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/27#issuecomment-868035006>.

- `gene_symbol`
- `gene_id`
- `cancer_group`
- `cohort`
- `tpm_mean`
- `tpm_sd`
- `tpm_mean_cancer_group_wise_zscore`/`tpm_mean_gene_wise_zscore`
- `tpm_mean_cancer_group_wise_quantiles`
- `n_samples`

### Results

#### All cohort summary statistics tables

The following wide tables are generated using the methods described above. Rows are genes. Columns are `cancer_group`s, except that the first two columns are gene symbol and gene Ensembl ID. If one gene symbol matches to multiple Ensembl IDs, the value of the Ensembl ID column is a comma separated list of all Ensembl IDs, e.g. `ENSG00000206952.3,ENSG00000281910.1`. In `inpiut/ens_symbol.tsv`, `SNORA50A` is mapped to both `ENSG00000206952.3` and `ENSG00000281910.1`.

- `results/cancer_group_all_cohort_mean_tpm.tsv`
- `results/cancer_group_all_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_all_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_all_cohort_mean_tpm_gene_wise_z_scores.tsv`

The samples used in each `cancer_group` are listed in `results/cancer_group_all_cohort_sample_metadata.tsv`. The columns are 1) `cancer_group`, 2) the number of samples in the `cancer_group`, and 3) the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `cancer_group`.

#### Individual cohort summary statistics tables

The following wide tables are generated using the methods described above. Rows are genes. Columns are `cancer_group_cohort`s, except that the first two columns are gene symbol and gene Ensembl ID. A `cancer_group_cohort` is a string that concatenates a `cancer_group` and a `cohort` by `___`, e.g. `Meningioma___PBTA`, `Neuroblastoma___GMKF`, and `Diffuse midline glioma___PBTA`. If one gene symbol matches to multiple Ensembl IDs, the value of the Ensembl ID column is a comma separated list of all Ensembl IDs, e.g. `ENSG00000206952.3,ENSG00000281910.1`. In `inpiut/ens_symbol.tsv`, `SNORA50A` is mapped to both `ENSG00000206952.3` and `ENSG00000281910.1`.

- `results/cancer_group_individual_cohort_mean_tpm.tsv`
- `results/cancer_group_individual_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_z_scores.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_cancer_group_wise_quantiles.tsv`
- `results/cancer_group_individual_cohort_mean_tpm_gene_wise_z_scores.tsv`

The samples used in each `cancer_group_cohort` are listed in `results/cancer_group_individual_cohort_sample_metadata.tsv`. The columns are 1) `cancer_group_cohort`, 2) the number of samples in the `cancer_group_cohort`, and 3) the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `cancer_group_cohort`.

#### Long summary statistics tables

The following long tables are generated using the wide tables.

- `long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz`
- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz`

Each row of the long table is a tab-delimited record of the following columns, as suggested by @jharenza and @taylordm at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/27#issuecomment-868035006>.

- `gene_symbol`
- `gene_id`
- `cancer_group`
- `cohort`
- `tpm_mean`
- `tpm_sd`
- `tpm_mean_cancer_group_wise_zscore`/`tpm_mean_gene_wise_zscore`
- `tpm_mean_cancer_group_wise_quantiles`
- `n_samples`

If the record is generated using all cohorts, the `cohort` column takes the value of `AllCohorts`.

Different from the wide tables, the `gene_id` does not contain comma-separated ENSG IDs. If one gene symbol matches to multiple Ensembl IDs, each Ensembl gene ID will become one row in the long table. For example:

- Wide table:

```text
# cancer_group_individual_cohort_mean_tpm.tsv
gene_symbol    gene_id    Adamantimomatous craniopharyngioma___PBTA    Atypical meningioma___PBTA    Atypical Teratoid Rhabdoid Tumor___PBTA
CDR1    ENSG00000184258.6,ENSG00000281508.1    111.1765    31.932000000000002    332.2896153846154
```

- Long table:

```text
# long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz
CDR1    ENSG00000184258.6    Adamantimomatous craniopharyngioma    PBTA    111.1765    181.70448082054318    0.11408800038186823    Highest expressed 25%    20
CDR1    ENSG00000281508.1    Adamantimomatous craniopharyngioma    PBTA    111.1765    181.70448082054318    0.11408800038186823    Highest expressed 25%    20
CDR1    ENSG00000184258.6    Atypical meningioma    PBTA    31.932000000000002    50.10366423725914    0.016467672549703015    Highest expressed 25%    5
CDR1    ENSG00000281508.1    Atypical meningioma    PBTA    31.932000000000002    50.10366423725914    0.016467672549703015    Highest expressed 25%    5
CDR1    ENSG00000184258.6    Atypical Teratoid Rhabdoid Tumor    PBTA    332.2896153846154    1075.378786603049    0.37106139758788453    Highest expressed 25%    26
CDR1    ENSG00000281508.1    Atypical Teratoid Rhabdoid Tumor    PBTA    332.2896153846154    1075.378786603049    0.37106139758788453    Highest expressed 25%    26
```

### Usage

1. Change working directory to local `OpenPBTA-analysis`.
2. Download data using `bash download-data.sh`. Make sure `data/gene-expression-rsem-tpm-collapsed.rds` and `data/histologies.tsv` are downloaded.
3. Run this analysis module in the continuous integration (CI) docker image using `./scripts/run_in_ci.sh bash analyses/rna-seq-expression-summary-stats/run-rna-seq-expression-summary-stats.sh`.

### Module structure

```text
.
├── 01-tpm-summary-stats.R
├── README.md
├── input
│   └── ens_symbol.tsv
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
│   ├── long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz
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
- `input/ens_symbol.tsv`: Gene symbol to ENSG ID conversion table. Shared by @kgaonkar6.

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

- `long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv`
- `long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv`
