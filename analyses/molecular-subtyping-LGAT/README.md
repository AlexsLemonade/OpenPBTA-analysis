# Molecular Subtyping of Low Grade Gliomaa

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)), Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

In this analysis we subtype LGAT samples according to the presence/absence of _BRAF_ fusions and BRAF V600E point mutations. 

### Preprocessing

The files in the `lgat-subset` were generated via `00-subset-files-for-LGAT.R` using files from the `release-v15-20200228`. If running locally using [the instructions below](#run-script), new `lgat-subset` files will be generated using the files symlinked in `data/`.

### Inputs from data download

* `pbta-histologies.tsv`: is used to subset samples where `short_histology == LGAT`
* `pbta-snv-consensus-mutation.maf.tsv.gz`: is used to find samples where `Hugo_Symbol == "BRAF" & HGVSp_Short == "p.V600E"`
* `pbta-fusion-putative-oncogenic.tsv`: is used to find samples with _BRAF_ fusions

### Run script

```sh
bash run_subtyping.sh
```

#### Order of scripts in analysis

`00-subset-files-for-LGAT.R`: generates subset of wgs LGAT sample annotated with `BRAF_V600E` column denoting presence (`Yes`) or absence (`No`) of BRAF V600E mutations.

`01-make-lgat-final-table.Rmd`: generates final table for LGAT subtyping from metadata, fusion, and consensus mutation data release files.
The logic for subtyping LGAT samples (defined by the combination of `sample_id` and `Kids_First_Participant_ID`) is as follows:

* If there is no fusion data and no mutation data available, a sample is labeled `LGG, To be classified`
* If no _BRAF_ fusion is detected and there is no mutation data available, a sample is labeled `LGG, To be classified`
* If there is no fusion data available or _BRAF_ fusion detected and a BRAF V600E mutation is present, a sample is labeled `LGG, BRAF V600E`
* If there is no fusion data available or _BRAF_ fusion detected and a BRAF V600E mutation is absent, a sample is labeled `LGG, BRAF wildtype`
* If there is a _BRAF_ fusion detected but no V600E mutation or no mutation data available, a sample is labeled `LGG, BRAF fusion`
* If there is a _BRAF_ fusion detected and a BRAF V600E mutation is present, a sample is labeled `LGG, BRAF fusion/V600E`
