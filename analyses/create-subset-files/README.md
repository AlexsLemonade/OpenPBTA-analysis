## Steps for creating subset files for CI

1. Update to the most recent release of the data by running `bash download-data.sh` in the root directory of the repository.
2. Run the `01-create_subset_files` notebook from the root directory of the repository:

```
Rscript -e "rmarkdown::render('analyses/create-subset-files/01-create_subset_files.Rmd', 
                              clean = TRUE)"
```
3. Ensure that only the output of `01-create_subset_files` is in `data/testing` (e.g., remove previous `md5sum.txt` files and `release-notes.md`. Navigate to `data/testing` and create the new `md5sum.txt` files:

```
cd data/testing
md5sum * > md5sum.txt
```
The files in `data/testing` are now ready to be uploaded to S3.
