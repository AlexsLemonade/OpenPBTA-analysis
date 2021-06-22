## Calculate TPM summary statistics within each cancer group and cohort

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

Within each cancer group and cohort, calculate TPM means, standard deviations, z-scores, and ranks.

### Methods

For each `cancer_group`, select one of the following two sets of samples:

- Samples from all `cohort`s, e.g. PBTA and GMKF.
- Samples from each individual `cohort`.

If >= 5 samples are selected, generate the following summary statistics:

- TPM means of each gene across all selected samples, and denote this vector as `mean_TPM_vector`.
- TPM standard deviations of each gene across all selected samples.
- z-scores of each gene across all genes, computed as `z_score_vector = (mean_TPM_vector - mean(mean_TPM_vector)) / sd(mean_TPM_vector)`.
- Ranks of mean TPM of each gene across all selected samples, which takes one of the following four values: `Highest expressed 25%`, `Expression between upper quartile and median`, `Expression between median and lower quartile`, and `Lowest expressed 25%`. If multiple genes have the same mean TPM value, their tied rank is the lowest rank, in order to be conservative on the description of their expression levels.

Combine each type of the summary statistics vectors into a table, with rows as genes, and columns as `cancer_group_cohort`.

### Results

#### All cohort summary statistics tables

The following tables are generated using the methods described above. Rows are genes. Columns are `cancer_group`s, except that the first two columns are gene symbol and gene Ensembl ID.

- `results/cancer_group_all_cohort_mean_tpm.tsv`
- `results/cancer_group_all_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_z_scores.tsv`
- `results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_quantiles.tsv`

The samples used in each `cancer_group` are listed in `results/cancer_group_all_cohort_sample_metadata.tsv`. The columns are 1) `cancer_group`, 2) the number of samples in the `cancer_group`, and 3) the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `cancer_group`.

#### Individual cohort summary statistics tables

The following tables are generated using the methods described above. Rows are genes. Columns are `cancer_group_cohort`s, except that the first two columns are gene symbol and gene Ensembl ID. A `cancer_group_cohort` is a string that concatenates a `cancer_group` and a `cohort` by `___`, e.g. `Meningioma___PBTA`, `Neuroblastoma___GMKF`, and `Diffuse midline glioma___PBTA`.

- `results/cancer_group_individual_cohort_mean_tpm.tsv`
- `results/cancer_group_individual_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_individual_cohort_cancer_group_wise_mean_tpm_z_scores.tsv`
- `results/cancer_group_individual_cohort_cancer_group_wise_mean_tpm_quantiles.tsv`

The samples used in each `cancer_group_cohort` are listed in `results/cancer_group_individual_cohort_sample_metadata.tsv`. The columns are 1) `cancer_group_cohort`, 2) the number of samples in the `cancer_group_cohort`, and 3) the comma separated list of `Kids_First_Biospecimen_ID`s of the samples in the `cancer_group_cohort`.

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
│   ├── cancer_group_all_cohort_cancer_group_wise_mean_tpm_quantiles.tsv
│   ├── cancer_group_all_cohort_cancer_group_wise_mean_tpm_z_scores.tsv
│   ├── cancer_group_all_cohort_mean_tpm.tsv
│   ├── cancer_group_all_cohort_sample_metadata.tsv
│   ├── cancer_group_all_cohort_standard_deviation_tpm.tsv
│   ├── cancer_group_individual_cohort_cancer_group_wise_mean_tpm_quantiles.tsv
│   ├── cancer_group_individual_cohort_cancer_group_wise_mean_tpm_z_scores.tsv
│   ├── cancer_group_individual_cohort_mean_tpm.tsv
│   ├── cancer_group_individual_cohort_sample_metadata.tsv
│   └── cancer_group_individual_cohort_standard_deviation_tpm.tsv
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
- `results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_z_scores.tsv`
- `results/cancer_group_all_cohort_cancer_group_wise_mean_tpm_quantiles.tsv`

- `results/cancer_group_individual_cohort_sample_metadata.tsv`
- `results/cancer_group_individual_cohort_mean_tpm.tsv`
- `results/cancer_group_individual_cohort_standard_deviation_tpm.tsv`
- `results/cancer_group_individual_cohort_cancer_group_wise_mean_tpm_z_scores.tsv`
- `results/cancer_group_individual_cohort_cancer_group_wise_mean_tpm_quantiles.tsv`
