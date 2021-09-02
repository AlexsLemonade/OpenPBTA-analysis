## Annotate CNV table with mutation frequencies

Adapted from [snv-frequencies](https://github.com/logstar/OpenPedCan-analysis/tree/snv-freq/analyses/snv-frequencies)
**Module author:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

Adapted from [fusion-frequencies](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/kgaonkar6/fusion_freq/analyses/fusion-frequencies)
**Module author:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6))

Adapted by Eric Wafula ([@ewafula](https://github.com/ewafula)) 

### Purpose
Uses `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz` and `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz` consensus CNV calls and variant types (`amplification`, `deep deletion`, `gain`, `loss`, and `neutral`) to determine `Ensembl` gene-level mutation frequencies for each cancer type in an overall cohort dateset and in the independent primary/relapse cohort subsets of the data.

#### Additional annotation
Additional disease and gene annotations include `gene full names` `RMTL designations`, `OncoKB categories`, and `EFO and MONDO identifiers` integrated to the CNV frequencies table using the [long-format-table-utils analysis module)](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/long-format-table-utils).

Each `cancer_group` and `cohort`(s) combination is considered a `cancer_group_cohort`. `cancer_group_cohort` with `n_samples` >= 5, compute `Frequency_in_overall_dataset`, `Frequency_in_primary_tumors`, and `Frequency_in_relapse_tumors` as following:

- `Frequency_in_overall_dataset`:
  - For each unique variant, count the number of patients (identified by `Kids_First_Participant_ID`) that have the variant, and call this number `Total_alterations`.
  - Count the total number of patients in the `cancer_group_cohort`, and call this number `Patients_in_dataset`.
  - `Frequency_in_overall_dataset = Total_alterations / Patients_in_dataset`.

- `Frequency_in_primary_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `independent-specimens.wgs.primary.eachcohort.tsv`, and call this number `Total_primary_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `independent-specimens.wgswxspanel.primary.eachcohort.tsv` for each cohort and `independent-specimens.wgswxspanel.primary.tsv` for all cohorts, and call this number `Primary_tumors_in_dataset`.
  - `Frequency_in_primary_tumors = Total_primary_tumors_alterated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique variant, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that are in the `independent-specimens.wgs.relapse.eachcohort.tsv`, and call this number `Total_relapse_tumors_alterated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the `independent-specimens.wgswxspanel.relapse.eachcohort.tsv` for each cohort and `independent-specimens.wgswxspanel.relapse.tsv` for all cohorts, and call this number `Relapse_tumors_in_dataset`.
  - `Frequency_in_relapse_tumors = Total_relapse_tumors_alterated / Relapse_tumors_in_dataset`.

Format the CNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/66>.

Merge the CNV frequency tables of all `cancer_group_cohort`s.

### Results

Results are generated using PediatricOpenTargets/OpenPedCan-analysis data release v7.

The merged CNV frequency table of all `cancer_group_cohort`s is output in TSV and JSONL formats.

- `gene-level-cnv-consensus-annotated-mut-freq.tsv.gz`
- `gene-level-cnv-consensus-annotated-mut-freq.jsonl.gz`

### Analysis scripts

### `run-cnv-frequencies-analysis.sh`
This is a bash script wrapper for setting input file paths for the main anlysis script, `01-cnv-frequencies.py`. All file paths in set in this script are based on root directory of this Git repository . Therefore, the script should always be run from the root directory of OPenPedCan-analysis

Adapted from [snv-callers analysis module](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/blob/dev/analyses/snv-callers/run_caller_consensus_analysis.sh)
**Module author:** Candace Savonen ([@cansavvy](https://github.com/cansavvy))

Usage:
```bash
bash analysis/cnv-frequencies/run-cnv-frequencies-analysis.sh

```

### `01-cnv-frequencies.py`
Python functions to create copy number variation (CNV) cancer type and study gene-level frequencies for OPenPedCan analyses modules

Usage:
```bash
python3 analysis/cnv-frequencies/01-cnv-frequencies.py HISTOLOGY_FILE CNV_FILE  AC_PRIMARY_TUMORS AC_RELAPSE_TUMORS EC_PRIMARY_TUMORS EC_RELAPSE_TUMORS
```

Parameter Options:
```
positional arguments:
  HISTOLOGY_FILE  OPenPedCan histology file (histologies.tsv)
                  
  CNV_FILE        OPenPedCan CNV consensus file 
                  (consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz
                  and consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz)
                  
  AC_PRIMARY_TUMORS  OPenPedCan all cohorts independent primary tumor samples file 
                  (independent-specimens.wgswxspanel.primary.tsv)
                  
  AC_RELAPSE_TUMORS  OPenPedCan all cohorts independent relapse tumor samples file 
                  (independent-specimens.wgswxspanel.relapse.tsv)
                  
  EC_PRIMARY_TUMORS  OPenPedCan each cohort independent primary tumor samples file 
                  (independent-specimens.wgswxspanel.primary.eachcohort.tsv)
                  
  EC_RELAPSE_TUMORS  OPenPedCan each cohort independent relapse tumor samples file 
                  (independent-specimens.wgswxspanel.relapse.eachcohort.tsv)

optional arguments:
  -h, --help      show this help message and exit
  -v, --version   Print the current 01-cnv-frequencies.py version and exit
```

Input:
- `data/histologies.tsv`
- `data/consensus_wgs_plus_cnvkit_wxs.tsv.gz`
- `data/independent-specimens.wgswxspanel.primary.tsv`
- `data/independent-specimens.wgswxspanel.relapse.tsv`
- `data/independent-specimens.wgswxspanel.primary.eachcohort.tsv`
- `data/independent-specimens.wgswxspanel.relapse.eachcohort.tsv`

Output:
- `analysis/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.tsv.gz`
- `analysis/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.jsonl.gz`

