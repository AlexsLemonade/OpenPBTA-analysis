## Annotate SNV table with mutation frequencies

**Module author:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

- [Annotate SNV table with mutation frequencies](#annotate-snv-table-with-mutation-frequencies)
  - [Purpose](#purpose)
  - [Methods](#methods)
    - [Subset overlapping tumor samples between `snv-consensus-plus-hotspots.maf.tsv.gz` and `histologies.tsv`](#subset-overlapping-tumor-samples-between-snv-consensus-plus-hotspotsmaftsvgz-and-histologiestsv)
    - [Subset non-synonymous variants](#subset-non-synonymous-variants)
    - [Generate variant-level and gene-level non-synonymous mutation frequencies](#generate-variant-level-and-gene-level-non-synonymous-mutation-frequencies)
    - [Add annotations](#add-annotations)
  - [Results](#results)
  - [Usage](#usage)
  - [Analysis scripts](#analysis-scripts)
    - [`01-snv-frequencies.R`](#01-snv-frequenciesr)
  - [Guide on how to update `01-snv-frequencies.R`](#guide-on-how-to-update-01-snv-frequenciesr)

### Purpose

- Annotate each non-synonymous variant in `snv-consensus-plus-hotspots.maf.tsv.gz` with mutation frequencies per `(cohort, cancer_group, primary/relapse)`.
- Annotate each gene in `snv-consensus-plus-hotspots.maf.tsv.gz` with non-synonymous variant mutation frequencies per `(cohort, cancer_group, primary/relapse)`.

Issues addressed:

- Variant-level mutation frequencies: <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>
- Gene-level mutation frequencies: <https://github.com/PediatricOpenTargets/ticket-tracker/issues/91>
- <https://github.com/PediatricOpenTargets/ticket-tracker/issues/8>. This issue is no longer compatible with the purpose of this module. This module intends to compute mutation frequencies for each variant, but this issue intents to compute the mutation frequencies for each gene. This issue is listed here for future reference.

### Methods

#### Subset overlapping tumor samples between `snv-consensus-plus-hotspots.maf.tsv.gz` and `histologies.tsv`

Subset `snv-consensus-plus-hotspots.maf.tsv.gz` to keep only samples with `sample_type == 'Tumor'` and non-NA `cancer_group` and `cohort` values in `histologies.tsv`.

Subset `histologies.tsv`, `../independent-samples/results/independent-specimens.wgswxspanel.primary.eachcohort.tsv`, and `../independent-samples/results/independent-specimens.wgswxspanel.relapse.eachcohort.tsv` to keep only samples that are in the `snv-consensus-plus-hotspots.maf.tsv.gz` subset.

#### Subset non-synonymous variants

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

#### Generate variant-level and gene-level non-synonymous mutation frequencies

Add `Gene_full_name` and `Protein_RefSeq_ID` columns to each variant with annotations obtained from [mygene.info](http://mygene.info/about).

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
  - For each unique variant/gene, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that have mutations at the variant/gene and are in the independent primary sample list, and call this number `Total_primary_tumors_mutated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the independent primary sample list, and call this number `Primary_tumors_in_dataset`.
  - `Frequency_in_primary_tumors = Total_primary_tumors_mutated / Primary_tumors_in_dataset`.

- `Frequency_in_relapse_tumors`:
  - For each unique variant/gene, count the number of samples (identified by `Kids_First_Biospecimen_ID`) that have mutations at the variant/gene and are in the independent relapse sample list, and call this number `Total_relapse_tumors_mutated`.
  - Count the total number of samples in the `cancer_group_cohort` that are also in the independent relapse sample list, and call this number `Relapse_tumors_in_dataset`.
  - `Frequency_in_relapse_tumors = Total_relapse_tumors_mutated / Relapse_tumors_in_dataset`.

Format the SNV mutation frequency table according to the latest spreadsheet that is attached in <https://github.com/PediatricOpenTargets/ticket-tracker/issues/64>.

Merge the SNV mutation frequency tables of all `cancer_group_cohort`s.

#### Add annotations

Add the following columns to variant-level and gene-level tables:

- EFO
- MONDO
- RMTL
- OncoKB cancer gene
- OncoKB oncogene/TSG (tumor suppresor gene)
- PedcBioPortal oncoprint plot URL
- PedcBioPortal lollipop plot URL columns

The EFO, MONDO, and RMTL information is obtained from PediatricOpenTargets/OpenPedCan-analysis data release.

The OncoKB cancer gene and oncogene/TSG is listed in `input/oncokb_cancer_gene_list.tsv`, which is downloaded from <https://www.oncokb.org/cancerGenes>. The last update of the table is on 06/16/2021. To update the table, re-download the updated table from <https://www.oncokb.org/cancerGenes>.

The PedcBioPortal `case_set_id`s in the URLs are obtained from [the `sample-lists` PedcBioPortal web API](https://pedcbioportal.kidsfirstdrc.org/api/swagger-ui.html#/Sample_Lists), with the following command:

```bash
curl -X GET "https://pedcbioportal.kidsfirstdrc.org/api/studies/ped_opentargets_2021/sample-lists" -H "accept: application/json" -H "Authorization: Bearer YOUR-API-ACCESS-TOKEN" > ped_opentargets_2021_pedcbio_case_set_ids.json
```

To update the `case_set_id`s, rerun the `curl` command with your own PedcBioPortal web API access token. The token can be requested and downloaded at <https://pedcbioportal.kidsfirstdrc.org/webAPI>. More information about the access token is at <https://docs.cbioportal.org/2.2-authorization-and-authentication/authenticating-users-via-tokens#using-data-access-tokens>.

Add a `Gene_type` column to gene-level tables. The `Gene_type` information is obtained from `../fusion_filtering/references/genelistreference.txt`, and its sources are described at <https://github.com/d3b-center/annoFuse#prerequisites-for-cohort-level-analysis>.

### Results

Results are generated using PediatricOpenTargets/OpenPedCan-analysis data release v6.

The merged variant-level and gene-level SNV mutation frequency tables of all `cancer_group_cohort`s is output in TSV, JSON, and JSONL formats.

- `results/gene-level-snv-consensus-annotated-mut-freq.json.gz`
- `results/gene-level-snv-consensus-annotated-mut-freq.jsonl.gz`
- `results/gene-level-snv-consensus-annotated-mut-freq.tsv`

- `results/var-level-snv-consensus-annotated-mut-freq.json.gz`
- `results/var-level-snv-consensus-annotated-mut-freq.jsonl.gz`
- `results/var-level-snv-consensus-annotated-mut-freq.tsv`

**Note:** When reading `{var,gene}-level-snv-consensus-annotated-mut-freq.tsv` using `readr::read_tsv`, you can specify `na = character(0)` to avoid converting blank string entries to `NA`s, as suggested by @jharenza at <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/45#pullrequestreview-697753423>. This can be useful when converting `{var,gene}-level-snv-consensus-annotated-mut-freq.tsv` to JSON, because the PediatricOpenTargets/OpenPedCan-analysis project favors empty string (`''`) as missing value over `NA`/`NaN`/`NULL`.

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
- `../../data/efo-mondo-map.tsv`
- `../../data/ensg-hugo-rmtl-v1-mapping.tsv`
- `../independent-samples/results/independent-specimens.wgs.primary.tsv`
- `../independent-samples/results/independent-specimens.wgs.relapse.tsv`
- `../fusion_filtering/references/genelistreference.txt`
- `input/oncokb_cancer_gene_list.tsv`
- `input/ped_opentargets_2021_pedcbio_case_set_ids.json`

Output:

- `results/gene-level-snv-consensus-annotated-mut-freq.json`
- `results/gene-level-snv-consensus-annotated-mut-freq.jsonl`
- `results/gene-level-snv-consensus-annotated-mut-freq.tsv`
- `results/var-level-snv-consensus-annotated-mut-freq.json`
- `results/var-level-snv-consensus-annotated-mut-freq.jsonl`
- `results/var-level-snv-consensus-annotated-mut-freq.tsv`

### Guide on how to update `01-snv-frequencies.R`

Read the script. `01-snv-frequencies.R` contains the following sections:

- Function definitions. This section contains the definitions of all functions used in this script. These functions directly use column names of input data, so any future updates need to check whether the column names are changed in new data releases. This script is originally developed using v6 data release.
- Create output dir. This section creates output directory.
- Read data. This section reads all input data.
- Subset tumor samples and used columns in MAF table. This section subsets MAF table.
- Subset independent samples in histology table. This section subsets histology table and creates primary and relapse independent sample histology tables.
- Add additional annotations. This section adds Gene_full_name and Protein_RefSeq_ID to the MAF table. Note that annotations will be added in the upcoming `analyses/long-format-table-utils/annotator/annotator-api.R`, which is under review at PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/56>.
- Compute mutation frequencies. This section:
  - computes mutation frequencies for 1) each cancer_group and each cohort and 2) each cancer_group and all cohorts
  - add each-cohot and all-cohort PedcBio URLs
  - creates gene-level and variant-level mutation frqeuency tables
- Add annotations to the output table. This section adds various annotations to the mutation frequency tables. Note that annotations will be added in the upcoming `analyses/long-format-table-utils/annotator/annotator-api.R`, which is under review at PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/56>.
- Output tsv and JSON. This section outputs mutation frequency tables in tsv and JSON formats.

Common updates:

- Add new gene annotations:
  - Load gene HUGO symbol (/gene Ensembl ENSG ID) annotation table in the Read data section. Make sure the HUGO symbol (/gene Ensembl ENSG ID) column has no NA or duplicate.
  - Add annotations to the mutation frequency tables in the Add annotations to the output table section. Replace `NA`s with `''`.
  - Reorder new columns before output.
  - Note that annotations will be added in the upcoming `analyses/long-format-table-utils/annotator/annotator-api.R`, which is under review at PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/56>.
- Add new variant annotations:
  - Load variant_ID annotation table in the Read data section. Make sure the variant_ID column has no NA or duplicate.
  - Add annotations to the MAF table in the Add additional annotations section. Replace `NA`s with `''`.
  - Retain the new annotation columns in the `get_cg_ch_var_level_mut_freq_tbl()` and `get_cg_ch_gene_level_mut_freq_tbl()` functions.
  - Reorder new columns before output.
  - Note that annotations will be added in the upcoming `analyses/long-format-table-utils/annotator/annotator-api.R`, which is under review at PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/56>.
- Change the notation for all_cohorts:
  - Change `'all_cohorts'` to other values in the `get_cohort_set_value()` function.
- Change the format of PedcBio URLs:
  - Change the `get_pcb_pot_csi()` and `get_pcb_pot_ploAt_url()` functions.
  - If `case_set_id`s are no longer used for linking PedcBio, the whole worflow needs to be redesigned.
- Add new PedcBio URLs:
  - Change the `add_cg_ch_pedcbio_pedot_plot_urls()` function.
  - Retain the new URL columns in the `get_cg_ch_var_level_mut_freq_tbl()` and `get_cg_ch_gene_level_mut_freq_tbl()` functions.
- Add new mutation frequency columns:
  - Change the `get_opr_mut_freq_tbl()` function.
  - Retain the new mutation frequency columns in the `get_cg_ch_var_level_mut_freq_tbl()` and `get_cg_ch_gene_level_mut_freq_tbl()` functions.

The `stopifnot()` statements are assertions for the input data, so the output would be expected. If any assertion fails, additional code needs to be added mainly to remove NAs and handle duplicates in the data, so that the assertion pases.

The unit testing is implemented using the [`testthat`](https://testthat.r-lib.org/index.html) package version 2.1.1, as suggested by @jharenza and @NHJohnson in the reviews of PR <https://github.com/PediatricOpenTargets/OpenPedCan-analysis/pull/55>.

To run all unit tests, run `bash run-tests.sh` in the Docker image/container from any working directory. Following is an example run.

```text
$ bash run-tests.sh
✔ |  OK F W S | Context
✔ |  19       | tests/test_collapse_rp_lists.R

══ Results ═════════════════════════════════════════════════════════════════════════════════════════════════════
Duration: 0.1 s

OK:       19
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
