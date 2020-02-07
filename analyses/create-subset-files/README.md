## Steps for creating subset files for CI

1. Update to the most recent release of the data by running `bash download-data.sh` in the root directory of the repository.
2. Run the shell script to generate subset files (from the root directory of the repository):

```
RELEASE=<RELEASE> ./analyses/create-subset-files/create_subset_files.sh
```
Alternatively, change the default value for `RELEASE` in the shell script.

3. The files in `data/testing/<RELEASE>` are now ready to be uploaded to S3.

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
