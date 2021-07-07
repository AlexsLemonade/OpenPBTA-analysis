## Annotate SNV table with mutation frequencies

**Module author:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

- Annotate each non-synonymous variant in `snv-consensus-plus-hotspots.maf.tsv.gz` with mutation frequencies per `(cohort, cancer_group, primary/relapse)`.
- Annotate each gene in `snv-consensus-plus-hotspots.maf.tsv.gz` with non-synonymous variant mutation frequencies per `(cohort, cancer_group, primary/relapse)`.

Issues addressed:

- Variant-level mutation frequencies: <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>
- Gene-level mutation frequencies: <https://github.com/PediatricOpenTargets/ticket-tracker/issues/91>
- <https://github.com/PediatricOpenTargets/ticket-tracker/issues/8>. This issue is no longer compatible with the purpose of this module. This module intends to compute mutation frequencies for each variant, but this issue intents to compute the mutation frequencies for each gene. This issue is listed here for future reference.

### Methods

Subset `snv-consensus-plus-hotspots.maf.tsv.gz` to keep only samples with `sample_type == 'Tumor'`.

Subset `snv-consensus-plus-hotspots.maf.tsv.gz` to keep only non-synonymous variants with the following code.

```R
Variant_Classification %in% c('Frame_Shift_Del',
                              'Frame_Shift_Ins',
                              'Splice_Site',
                              'Nonsense_Mutation',
                              'Nonstop_Mutation',
                              'In_Frame_Del',
                              'In_Frame_Ins',
                              'Missense_Mutation',
                              'Fusion',
                              'Multi_Hit',
                              'Multi_Hit_Fusion',
                              'Hom_Deletion',
                              'Hem_Deletion',
                              'Amp',
                              'Del',
                              'Translation_Start_Site')
```

Add `Gene_full_name` and `Protein_RefSeq_ID` columns to each variant with annotations obtained from [mygene.info](http://mygene.info/about).

Generate variant-level and gene-level non-synonymous mutation frequencies using the following procedures.

For vairant-level analysis, create a `Variant_ID` for each variant by concatenating `Chromosome`, `Start_Position`, `Reference_Allele`, and `Tumor_Seq_Allele2` with `'_'`.

For each `cancer_group`, get each cohort and all cohorts. Call each `cancer_group` and `cohort`(s) combination as a `cancer_group_cohort`. For example,

| cancer_group  | cohort           | n_samples |
|---------------|------------------|-----------|
| Neuroblastoma | CBTN             | 2         |
| Neuroblastoma | GMKF             | 541       |
| Neuroblastoma | TARGET           | 889       |
| Neuroblastoma | CBTN&GMKF&TARGET | 1432      |

For each `cancer_group_cohort` with `n_samples` >= 5, compute `Frequency_in_overall_dataset`, `Frequency_in_primary_tumors`, and `Frequency_in_relapse_tumors` as following:

- `Frequency_in_overall_dataset`:
  - For each unique variant/gene, count the number of patients (identified by `Kids_First_Participant_ID`) that have mutations at the variant/gene, and call this number `Total_mutations`.
  - Count the total number of patients in the `cancer_group_cohort`, and call this number `Patients_in_dataset`.
  - `Frequency_in_overall_dataset = Total_mutations / Patients_in_dataset`.

- `Frequency_in_primary_tumors`:
  - For each unique variant/gene, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that have mutations at the variant/gene and are in the `../independent-samples/results/independent-specimens.wgs.primary.tsv`, and call this number `Total_primary_tumors_mutated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.wgs.primary.tsv`, and call this number `Primary_tumors_in_dataset`.
  - `Frequency_in_primary_tumors = Total_primary_tumors_mutated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique variant/gene, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that have mutations at the variant/gene and are in the `../independent-samples/results/independent-specimens.wgs.relapse.tsv`, and call this number `Total_relapse_tumors_mutated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.wgs.relapse.tsv`, and call this number `Relapse_tumors_in_dataset`.
  - `Frequency_in_relapse_tumors = Total_relapse_tumors_mutated / Relapse_tumors_in_dataset`.

Format the SNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>.

Merge the SNV mutation frequency tables of all `cancer_group_cohort`s.

### Results

Results are generated using PediatricOpenTargets/OpenPedCan-analysis data release v6.

The merged SNV mutation frequency table of all `cancer_group_cohort`s is output in TSV and JSON formats using `readr::write_tsv()` and `jsonlite::write_json()` respectively.

- `results/snv-consensus-annotated-mut-freq.tsv`
- `snv-consensus-annotated-mut-freq.json.gz`

**Note:** When reading `snv-consensus-annotated-mut-freq.tsv` using `readr::read_tsv`, you can specify `na = character(0)` to avoid converting blank string entries to `NA`s, as suggested by @jharenza at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/45#pullrequestreview-697753423>. This can be useful when converting `snv-consensus-annotated-mut-freq.tsv` to JSON, because the PediatricOpenTargets/OpenPedCan-analysis project favors empty string (`''`) as missing value over `NA`/`NaN`/`NULL`.

### Usage

1. Change working directory to local `OpenPBTA-analysis`.
2. Download data using `bash download-data.sh`.
3. Run this analysis module in the continuous integration (CI) docker image using `./scripts/run_in_ci.sh bash analyses/snv-frequencies/run-snv-frequencies.sh`.

### Analysis scripts

#### `01-snv-frequencies.R`

This script annotates each non-synonymous variant in `snv-consensus-plus-hotspots.maf.tsv.gz` with mutation frequencies per `(cohort, cancer_group, primary/relapse)`.

Usage:

```bash
Rscript --vanilla '01-snv-frequencies.R'
```

Input:

- `../../data/histologies.tsv`
- `../../data/snv-consensus-plus-hotspots.maf.tsv.gz`
- `../independent-samples/results/independent-specimens.wgs.primary.tsv`
- `../independent-samples/results/independent-specimens.wgs.relapse.tsv`
- `../../data/efo-mondo-map.tsv`
- `../../data/ensg-hugo-rmtl-v1-mapping.tsv`
- `../fusion_filtering/references/genelistreference.txt`

Output:

- `results/snv-consensus-annotated-mut-freq.tsv`
- `results/snv-consensus-annotated-mut-freq.json`
