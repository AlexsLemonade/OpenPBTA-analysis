# Molecular Subtyping of Low Grade Gliomaa

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

In this analysis we subtype LGAT samples according to the presence/absence of _BRAF_ fusions and BRAF V600E point mutations. 

### Inclusion/exclusion criteria

Samples are _included_ for subtyping if we detect the following strings in the `pathology_diagnosis` field of `pbta-histologies.tsv`:

```
low-grade glioma/astrocytoma
ganglioglioma
```

Samples are _excluded_ if we detect the following strings in the `pathology_diagnosis` field of `pbta-histologies.tsv`:

```
dysembryoplastic neuroepithelial tumor
```

These strings are stored in `lgat-subset/lgat_subtyping_path_dx_strings.json`, which is the output of the `00-LGAT-select-pathology-dx` notebook. 

Prior to `release-v17-20200908`, we used `short_histology == "LGAT"` as our inclusion criteria.


### Preprocessing

The files in the `lgat-subset` were generated via `01-subset-files-for-LGAT.R` using files from the `release-v17-20200908`. If running locally using [the instructions below](#run-script), new `lgat-subset` files will be generated using the files symlinked in `data/`.

### Inputs from data download

* `pbta-histologies.tsv`: is used to subset samples according to [the criteria above](#inclusion-exclusion-criteria)
* `pbta-snv-consensus-mutation.maf.tsv.gz`: is used to find samples where `Hugo_Symbol == "BRAF" & HGVSp_Short == "p.V600E"`
* `pbta-fusion-putative-oncogenic.tsv`: is used to find samples with _BRAF_ fusions

### Run script

```sh
bash run_subtyping.sh
```

This does not run the `00-LGAT-select-pathology-dx` notebook, as that is intended to be run once and tied to a specific release (`release-v17-20200908`).

#### Order of scripts in subtyping

`01-subset-files-for-LGAT.R`: generates subset of wgs LGAT sample annotated with `BRAF_V600E` column denoting presence (`Yes`) or absence (`No`) of BRAF V600E mutations.

`02-make-lgat-final-table.Rmd`: generates final table for LGAT subtyping from metadata, fusion, and consensus mutation data release files.
The logic for subtyping LGAT samples (defined by the combination of `sample_id` and `Kids_First_Participant_ID`) is as follows:

* If there is no fusion data and no mutation data available, a sample is labeled `LGG, To be classified`
* If no _BRAF_ fusion is detected and there is no mutation data available, a sample is labeled `LGG, To be classified`
* If there is no fusion data available or _BRAF_ fusion detected and a BRAF V600E mutation is present, a sample is labeled `LGG, BRAF V600E`
* If there is no fusion data available or _BRAF_ fusion detected and a BRAF V600E mutation is absent, a sample is labeled `LGG, BRAF wildtype`
* If there is a _BRAF_ fusion detected but no V600E mutation or no mutation data available, a sample is labeled `LGG, BRAF fusion`
* If there is a _BRAF_ fusion detected and a BRAF V600E mutation is present, a sample is labeled `LGG, BRAF fusion/V600E`
