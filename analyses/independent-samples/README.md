# Independent Samples

## Module authors
Komal Rathi, Run Jin, Yuanchao Zhang

## Module structure

* `01-generate-independent-specimens-wgs-only.R`: Generate tables of WGS-only independent specimens where no two specimens are chosen from the same individual.
* `01-generate-independent-specimens-wgs-preferred.R`: Generate tables of WGS-preferred independent specimens where no two specimens are chosen from the same individual. 
* `01-generate-independent-specimens-wxs-preferred.R`: Generate tables of WXS-preferred independent specimens where no two specimens are chosen from the same individual.
* `02-generate-independent-rnaseq.R`: Generate tables of independent rna-seq specimens.
* `03-qc-independent-samples.Rmd`: Markdown to tabulate number of biospecimen ids for same participant ids from each output file.

```
.
├── 00-repeated-samples.Rmd
├── 00-repeated-samples.nb.html
├── 01-generate-independent-specimens-wgs-only.R
├── 01-generate-independent-specimens-wgs-preferred.R
├── 01-generate-independent-specimens-wxs-preferred.R
├── 02-generate-independent-rnaseq.R
├── 03-qc-independent-samples.Rmd
├── 03-qc-independent-samples.nb.html
├── README.md
├── results
│   ├── independent-specimens.rnaseq.primary-plus.eachcohort.tsv
│   ├── independent-specimens.rnaseq.primary-plus.tsv
│   ├── independent-specimens.rnaseq.primary.eachcohort.tsv
│   ├── independent-specimens.rnaseq.primary.tsv
│   ├── independent-specimens.rnaseq.relapse.eachcohort.tsv
│   ├── independent-specimens.rnaseq.relapse.tsv
│   ├── independent-specimens.wgs.primary-plus.eachcohort.tsv
│   ├── independent-specimens.wgs.primary-plus.tsv
│   ├── independent-specimens.wgs.primary.eachcohort.tsv
│   ├── independent-specimens.wgs.primary.tsv
│   ├── independent-specimens.wgs.relapse.eachcohort.tsv
│   ├── independent-specimens.wgs.relapse.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.tsv
│   ├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.primary.eachcohort.tsv
│   ├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.primary.tsv
│   ├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.relapse.eachcohort.tsv
│   ├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
│   └── independent-specimens.wgswxspanel.relapse.tsv
├── run-independent-samples.sh
└── util
    ├── independent-samples.R
    └── independent_rna_samples.R
```

## Summary

Many analyses that involve mutation frequencies or co-occurence require that all samples be independent. However, the `OpenPedCan` data set includes many cases where multiple biospecimens were taken from a single participant. This analysis creates lists of samples such that there are no cases where more than one biospecimen is included from each participant.

As different analyses may require different sets of data, we generate a few different sets, stored in the `results` subdirectory. We run the analyses based on different `independent_level`, either `each-cohort` or `all-cohorts`. For `each-cohort` analysis, we generate a list of independent samples using `Kids_First_Participant_ID` unique within each `cohort + cancer_group` combination. In this case, same `Kids_First_Participant_ID` across different cohorts is treated as `independent`. For `all-cohorts` analysis, we generate a list of independent samples using `Kids_First_Participant_ID` regardless of `cohort` or `cancer_group` and only one occurence of `Kids_First_Participant_ID` is chosen from the entire dataset.

## Outputs

### WGS-only lists

These lists contain only WGS samples:

1. **All-cohorts specific lists**

* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.tsv`
* Relapse specimens with WGS:  
`independent-specimens.wgs.relapse.tsv`
* Primary and relapse specimens with WGS:  
`independent-specimens.wgs.primary-plus.tsv`

2. **Each-cohort specific lists**

* Primary specimens only with whole genome sequence (WGS):  
`independent-specimens.wgs.primary.eachcohort.tsv`
* Relapse specimens with WGS:  
`independent-specimens.wgs.relapse.eachcohort.tsv`
* Primary and relapse specimens with WGS:  
`independent-specimens.wgs.primary-plus.eachcohort.tsv`

### WGS-preferred lists

For WGS-preferred lists, when a `Kids_First_Participant_ID` is associated with multiple `experimental_strategy` values i.e. `WGS`, `WXS` or `Targeted Sequencing`, priority is given to a single randomly chosen `WGS` biospecimen first, followed by either a single randomly chosen `WXS` or `Targeted Sequencing` sample.

Here, we first subset the tumor samples to `WGS` samples only and generate `WGS-specific` lists. These lists only contain a single occurence of `Kids_First_Participant_ID` associated to the `experimental_strategy = WGS`. Next, we subset the tumor samples to `WXS and Targeted Sequencing` to generate `WXS/Panel specific` lists. These lists only contain a single occurence of `Kids_First_Participant_ID` associated to either `experimental_strategy = WXS` or `experimental_strategy = Targeted Sequencing`. Then we merge the two lists (keeping `WGS` list first and `WXS/Panel` as second) and take a `dplyr::distinct` to only get the first occurence of `Kids_First_Participant_ID`. Because we keep the `WGS` specific list first when calling `dplyr::distinct`, the `WGS` associated biospecimens will be preferred over other biospecimens for multiple occurrences of `Kids_First_Participant_ID`.

1. **All-cohorts specific lists**

* Primary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.primary.tsv`
* Relapse specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.relapse.tsv`
* Primary and relapse specimens with WGS or WXS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.tsv`

2. **Each-cohort specific lists**

* Primary specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.primary.eachcohort.tsv`
* Relapse specimens only with either WGS or whole exome sequence (WXS) or Panel:  
`independent-specimens.wgswxspanel.relapse.eachcohort.tsv`
* Primary and relapse specimens with WGS or WXS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv`

### WXS-preferred lists

Additionally, similar independent lists that consist of all DNA `experimental_strategy` were generated prioriting `WXS` specimens when both WXS and WGS samples were available - this independent list will be used for analyzing SNV datasets since WXS gives higher coverage for coding regions and hence, better chance of detecting SNV.

For WXS-preferred lists, when a `Kids_First_Participant_ID` is associated with multiple `experimental_strategy` values i.e. `WGS`, `WXS` or `Targeted Sequencing`, priority is given to a single randomly chosen `WXS` biospecimen first, followed by either a single randomly chosen `WGS` or `Targeted Sequencing` sample.

Here, we first subset the tumor samples to `WXS` samples only and generate `WXS-specific` lists. These lists only contain a single occurence of `Kids_First_Participant_ID` associated to the `experimental_strategy = WXS`. Next, we subset the tumor samples to `WGS and Targeted Sequencing` to generate `WGS/Panel specific` lists. These lists only contain a single occurence of `Kids_First_Participant_ID` associated to either `experimental_strategy = WGS` or `experimental_strategy = Targeted Sequencing`. Then we merge the two lists (keeping `WXS` list first and `WGS/Panel` as second) and take a `dplyr::distinct` to only get the first occurence of `Kids_First_Participant_ID`. Because we keep the `WXS` specific list first when calling `dplyr::distinct`, the `WXS` associated biospecimens will be preferred over other biospecimens for multiple occurrences of `Kids_First_Participant_ID`.

1. **All-cohorts specific lists**

* Primary specimens only with either whole exome sequence (WXS) or WGS or Panel:  
`independent-specimens.wgswxspanel.primary.prefer.wxs.tsv`
* Relapse specimens only with either whole exome sequence (WXS) or WGS or Panel:  
`independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv`
* Primary and relapse specimens with WXS or WGS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv`

2. **Each-cohort specific lists**

* Primary specimens only with either whole exome sequence (WXS) or WGS or Panel:  
`independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv`
* Relapse specimens only with either whole exome sequence (WXS) or WGS or Panel:  
`independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv`
* Primary and relapse specimens with WXS or WGS or Panel:  
`independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv`

### RNA-sequencing lists

Simiarly, for independent RNA samples, we also run with either `all-cohorts` or `each-cohort`.
When run with `each-cohort`, independent DNA samples ran with `each-cohort` was used as starting point (see code for details) and when run with `all-cohorts`, independent DNA samples ran with `all-cohorts` was used as starting point.

When multiple RNA-Seq samples exist per participant, the script matches the independent whole genome or whole exome sample_ids to gather matched RNA-Seq sample. If participant has only RNA-Seq sample then a primary (and relapse if applicable) sample is randomly selected per participant per cancer group per cohort. 

1. **All-cohorts specific lists**

* Primary RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary-plus.tsv`
* Relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary.tsv`
* Primary and Relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.relapse.tsv`

2. **Each-cohort specific lists**

* Primary RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary-plus.eachcohort.tsv`
* Relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.primary.eachcohort.tsv`
* Primary and Relapse RNA-Seq specimens matching WGS/WXS/Panel independent sample_ids plus only-RNA-Seq 
`independent-specimens.rnaseq.relapse.eachcohort.tsv`

## Generating sample lists

To generate the independent sample lists and associated analysis of redundancies in the overall data set, run the following script from the project root directory:

```sh
bash analyses/independent-samples/run-independent-samples.sh
```

## Additional info:
- When presented with more than one specimen from a given individual with a specific cancer group and cohort, the script selects the first occurence of the individual so as to include only one specimen, with preference for primary tumors and whole genome sequences where available.
- The input histology file is randomized using a seed before using as input in order to avoid any selection bias. Using a seed allows for reproducibility of the randomized histology file.
- There is also a preference for the earliest collected samples, but as this data is not currently available, that code is currently deleted.

