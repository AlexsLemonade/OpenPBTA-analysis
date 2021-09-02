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

For each `cancer_group` and `cohort`(s) combination is considered a `cancer_group_cohort`. `cancer_group_cohort` with `n_samples` >= 5,  `Frequency_in_overall_dataset`, `Frequency_in_primary_tumors`, and `Frequency_in_relapse_tumors` were computed. 
In this module, we generated two fusion frequency tables - one `putative-oncogene-fusion-freq.jsonl.gz` calculate the frequency of a specific fusion in primary tumors, relapse tumors, and overall dataset. 
In `putative-oncogene-fused-gene-freq.jsonl.gz` on the other hand, we calculate the frequency of a particular gene being fused in primary tumors, relapse tumors, and overall dataset.
As long as a gene is a partner in a particular fusion (from 1A, 1B, 2A or 2B), it is counted as fused in that particular sample. 

Table one - `putative-oncogene-fusion-freq.jsonl.gz`
- `Frequency_in_overall_dataset`:
  - For each unique fusion, count the number of patients (identified by `Kids_First_Participant_ID`) that have the fusion, and call this number `Total_alterations`.
- Count the total number of patients in the `cancer_group_cohort`, and call this number `Patients_in_dataset`.
- `Frequency_in_overall_dataset = Total_alterations / Patients_in_dataset`.

Please **NOTE** that `Total_alterations` is the number of patients that have the particular fusion - NOT the total number of fusions. 
If a specific fusion occured 3 times but they were all from the same patient, the `Total_alterations` would equal to 1 not 3.
Therefore, `Total_alterations` cannot be directly summed in most cases, because the value corresponds to the size of a set of patients, and union of the patient sets need to be taken to have a grand "sum" of Total_alterations. For example,
- Row 1 altered patients are {a, b, c}; Total_alterations = 3.
- Row 2 altered patients are {a, c, d}; Total_alterations = 3.
- Row 1 and 2 altered patients should be {a, b, c, d}; Total_alterations = 4.

- `Frequency_in_primary_tumors`:
  - For each unique fusion, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`, and call this number `Total_primary_tumors_alterated`.
- Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.rnaseq.primary.eachcohort.tsv`, and call this number `Primary_tumors_in_dataset`.
- `Frequency_in_primary_tumors = Total_primary_tumors_alterated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique fusion, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `../independent-samples/results/independent-specimens.rnaseq.relapase.eachcohort.tsv`, and call this number `Total_relapse_tumors_alterated`.
- Count the total number of samples in the `cancer_group_cohort` that are also in the `../independent-samples/results/independent-specimens.rnaseq.relapase.eachcohort.tsv`, and call this number `Relapse_tumors_in_dataset`.
- `Frequency_in_relapse_tumors = Total_relapse_tumors_alterated / Relapse_tumors_in_dataset`.

Please **NOTE** that just like `Total_alterations`, `Total_relapse_tumors_alterated` and `Total_primary_tumors_alterated` are also the number of patients that have the particular fusion - NOT the total number of fusions. 

Format the SNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>.

Merge the fusion frequency tables of all `cancer_group_cohort`s.
The combined table was is written out as `putative-oncogene-fusion-freq.jsonl.gz`.

Table two - `putative-oncogene-fused-gene-freq.jsonl.gz`
The overall calculation logic is identical to table one. 
The difference is instead of calculating the frequency for each unique fusion, we are calculating the frequency for each unique gene that is present within a fusion (from 1A, 1B, 2A or 2B in a particular fusion) in different tumor types.'

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
