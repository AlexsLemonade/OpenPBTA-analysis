## Annotate SNV table with mutation frequencies

Adapted from [snv-frequencies](https://github.com/logstar/OpenPedCan-analysis/tree/snv-freq/analyses/snv-frequencies)
**Module author:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

Adapted by Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)) 

### Purpose
Uses `fusion-putative-oncogenic.tsv` and a `Alt_ID` for each `FusionName` and `Fusion_Type` concatenated by "_" OR Fused Gene to count occurence in  each cancer_group in a given dataset, primary or relapse cohorts.

Each gene invloved in the fusion is annotated by a Gene_Position, for example:
 - In a genic fusion where both breakpoints are within gene body the Gene_Position `Gene1A` will be the 5' gene and Gene1B will be 3' gene
 - In an intergenic fusion where one or both breakpoint are outside the gene body, if the 5' breakpoint is a region between GeneX and GeneY  the Gene_Position of GeneX will be Gene1A, for GeneY will be Gene2A and their 3' breakpoint is within GeneE will be denoted as Gene_Position `Gene1B`.
 - In an intergenic fusion where one or both breakpoint are outside the gene body, if the 3' breakpoint is within GeneE the Gene_Position of GeneE will be Gene1A and the 5' breakpoint is within GeneY and GeneZ, Gene Y will be denoted as Gene_Position `Gene1B` and GeneZ will be denoted as `Gene2B`.   


#### Frequency annotation

Each `cancer_group` and `cohort`(s) combination is considered a `cancer_group_cohort`. `cancer_group_cohort` with `n_samples` >= 5, compute `Frequency_in_overall_dataset`, `Frequency_in_primary_tumors`, and `Frequency_in_relapse_tumors` as following:

- `Frequency_in_overall_dataset`:
  - For each unique variant, count the number of patients (identified by `Kids_First_Participant_ID`) that have the variant, and call this number `Total_alterations`.
  - Count the total number of patients in the `cancer_group_cohort`, and call this number `Patients_in_dataset`.
  - `Frequency_in_overall_dataset = Total_alterations / Patients_in_dataset`.

- `Frequency_in_primary_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`, and call this number `Total_primary_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`, and call this number `Primary_tumors_in_dataset`.
  - `Frequency_in_primary_tumors = Total_primary_tumors_alterated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.rnaseq.relapase.eachcohort.tsv`, and call this number `Total_relapse_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.rnaseq.relapase.eachcohort.tsv`, and call this number `Relapse_tumors_in_dataset`.
  - `Frequency_in_relapse_tumors = Total_relapse_tumors_alterated / Relapse_tumors_in_dataset`.

Format the SNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>.

Merge the fusion frequency tables of all `cancer_group_cohort`s.

#### Additional annotation

[long-format-table-utils](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/long-format-table-utils) provides the annotation for columns  EFO, MONDO, RMTL and Gene_full_name.

### Analysis scripts

### `01-fusion-frequencies.R`
This script annotates each FusionName with occurence of Fusion_Type OR Fused Gene and frequencies in each cancer_cohort in dataset, primary and relapse.


Usage:

```bash
Rscript --vanilla '01-fusion-frequencies.R'

```

Input:

- `../../data/histologies.tsv`
- `../fusion_filtering/results/fusion-putative-oncogenic.tsv`
- `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`
- `../independent-samples/results/independent-specimens.rnaseq.relapse.eachcohort.tsv`

```
results/
├── putative-oncogene-fused-gene-freq.json.gz
├── putative-oncogene-fused-gene-freq.jsonl.gz
├── putative-oncogene-fused-gene-freq.tsv.gz
├── putative-oncogene-fusion-freq.json.gz
├── putative-oncogene-fusion-freq.jsonl.gz
└── putative-oncogene-fusion-freq.tsv.gz
```
