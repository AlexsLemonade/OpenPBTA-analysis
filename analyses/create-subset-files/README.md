## Steps for creating subset files for CI

1. Update to the most recent release of the data by running `bash download-data.sh` in the root directory of the repository.
2. Run the shell script to generate subset files (from the root directory of the repository):

```
RELEASE=<RELEASE> ./analyses/create-subset-files/create_subset_files.sh
```
Alternatively, change the default value for `RELEASE` in the shell script.

3. The files in `data/testing/<RELEASE>` are now ready to be uploaded to S3.

### Overview of approach

The objective of this module is to create files that are used in continuous integration to save on the amount of time it takes to download files and run analyses as well as reducing the amount of RAM that is required.

We create these files by randomly selecting participants that are represented in each file. 
The number of matched participants (participants that will be represented in _all_ files in CI) is controlled by the `--num_matched` to `01-get_biospecimen_identifiers.R` (Rscript default is `25`) and the `NUM_MATCHED` environmental variable in `create_subset_files.sh` (the default is `15` and therefore this is the default for the module).
Non-matched samples are also added to each file (10% of `--num_matched`), which allows us to test the realistic scenario where a participant ID is in one file but not another.

Some files are copied over in their entirety (e.g., BED files).
See `create_subset_files.sh` for more information.

#### Special considerations

Certain analysis modules have required modifications to the subset file creation steps beyond randomly selecting participants.

* We have two RNA-seq datasets. The smaller of the two, `polya`, is overrepresented relative to the total proportion of the RNA-seq samples it contains.
* [`sex-prediction-from-RNASeq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/sex-prediction-from-RNASeq) required that we specified `Male` and `Female` samples (`reported_gender`) were both represented in the `stranded` RNA-seq data.
* `tp53_nf1_module` requires us to include a set of biospecimen IDs for samples that have _TP53_ and _NF1_ mutations and are present in the `stranded` RNA-seq dataset.
See the `00-enrich-positive-examples` notebook for more information.
* `fusion-summary` requires us to include fusions that involve _MN1_, _RELA_, or _EWSR1_.

### Local development

Some of the SNV files are quite large and are not amenable to subsetting with less than 128 GB of RAM.

To skip the two larger MAF files (Vardict, Mutect2), you can set the `--local` option of `01-get_biospecimen_identifiers.R` to `1`. 

To run the entire pipeline with skipping those files enabled, one can run (from the root directory of the repository):

```
RUN_LOCAL=1 ./analyses/create-subset-files/create_subset_files.sh
```

### Additional features for convenience

Running the following from the root directory of the repository

```
SKIP_SUBSETTING=1 ./analyses/create-subset-files/create_subset_files.sh
```

will skip the subsetting file steps that are implemented in R and only copy files that are included in full (e.g., `pbta-histologies.tsv`) and generate a new `md5sum.txt`.
This is intended to be used when the only files that need to be updated are those that are copied over without being reduced in size in anyway.

