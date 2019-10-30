## Steps for creating subset files for CI

1. Update to the most recent release of the data by running `bash download-data.sh` in the root directory of the repository.
2. Run the shell script to generate subset files (from the root directory of the repository):

```
RELEASE=<RELEASE> ./analyses/create-subset-files/create_subset_files.sh
```
Alternatively, change the default value for `RELEASE` in the shell script.

3. The files in `data/testing/<RELEASE>` are now ready to be uploaded to S3.
