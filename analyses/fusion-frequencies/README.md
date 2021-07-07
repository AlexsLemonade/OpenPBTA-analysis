## Annotate SNV table with mutation frequencies

Adapted from [snv-frequencies](https://github.com/logstar/OpenPedCan-analysis/tree/snv-freq/analyses/snv-frequencies)
**Module author:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

Adapted by Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)) 

### Purpose
Uses `fusion-putative-oncogenic.tsv` and a `Alt_ID` for each `FusionName` and `Fusion_Type` concatenated by "_" to count occurence in  each cancer_group in a given dataset, primary or relapse cohorts.


#### Additional annotation

Add `Gene_full_name` and `Protein_RefSeq_ID` columns to each variant with annotations obtained from [mygene.info](http://mygene.info/about).

Each `cancer_group` and `cohort`(s) combination is considered a `cancer_group_cohort`. `cancer_group_cohort` with `n_samples` >= 5, compute `Frequency_in_overall_dataset`, `Frequency_in_primary_tumors`, and `Frequency_in_relapse_tumors` as following:

- `Frequency_in_overall_dataset`:
  - For each unique variant, count the number of patients (identified by `Kids_First_Participant_ID`) that have the variant, and call this number `Total_alterations`.
  - Count the total number of patients in the `cancer_group_cohort`, and call this number `Patients_in_dataset`.
  - `Frequency_in_overall_dataset = Total_alterations / Patients_in_dataset`.

- `Frequency_in_primary_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.wgs.primary.tsv`, and call this number `Total_primary_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.wgs.primary.tsv`, and call this number `Primary_tumors_in_dataset`.
  - `Frequency_in_primary_tumors = Total_primary_tumors_alterated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.wgs.relapse.tsv`, and call this number `Total_relapse_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.wgs.relapse.tsv`, and call this number `Relapse_tumors_in_dataset`.
  - `Frequency_in_relapse_tumors = Total_relapse_tumors_alterated / Relapse_tumors_in_dataset`.

Format the SNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>.

Merge the SNV mutation frequency tables of all `cancer_group_cohort`s.

### Results

Results are generated using PediatricOpenTargets/OpenPedCan-analysis data release v5.

`results/snv-consensus-annotated-mut-freq.tsv` is the merged Fusion frequency table for all `cancer_group_cohort`s

### Analysis scripts

### `fusion-frequencies.R`
This script annotates each FusionName with occurence of Fusion_Type and frequencies in each cancer_cohort in dataset, primary and relapse.


Usage:

```bash
Rscript --vanilla '01-fusion-frequencies.R'

```

Input:

- `../../data/histologies.tsv`
- `../fusion_filtering/results/fusion-putative-oncogenic.tsv`
- `../independent-samples/results/independent-specimens.wgs.primary.tsv`
- `../independent-samples/results/independent-specimens.wgs.relapse.tsv`

Output:
- putative-oncogene-fusion-annotated-freq.json
- putative-oncogene-fusion-annotated-freq.tsv
