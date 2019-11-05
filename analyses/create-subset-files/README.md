## Creating subset files for CI

**The first step is to update to the most recent release of the data by running `bash download-data.sh` in the root directory of the repository.**

### Method 1: Generate all subset files, including new biospecimen ID lists

Use this method when new _data_ files have been added to the data release or samples have been added or hidden from files.
This can always be run, but it is time- and RAM-intensive.

To use this method, run the following from the root directory of the repository:

```
RELEASE=<RELEASE> ./analyses/create-subset-files/create_all_subset_files.sh
```

See `create_all_subset_files.sh` for additional arguments.
The files for CI will be in `data/testing/<RELEASE>`.

### Method 2: Generate selected subset files when older subset files are available

Use this method if you have a local copy of the subset files for an older release (i.e., you have the `data/testing/<older_release>`) and you need to generate subset files for a selection of files that were updated in the new release.

To use this method, run the `create_selected_subset_files.sh` script from the root directory.
Here is an example for how this would be run for the files changed between `release-v7-20191031` and `release-v8-20191104`:

```
OLD_RELEASE=release-v7-20191031 NEW_RELEASE=release-v8-20191104 SELECTED_FILES='pbta-cnv-controlfreec.tsv.gz,pbta-gene-counts-rsem-expected_count.polya.rds,pbta-gene-counts-rsem-expected_count.stranded.rds,pbta-gene-expression-rsem-fpkm.polya.rds,pbta-gene-expression-rsem-fpkm.stranded.rds,pbta-histologies.tsv' ./analyses/create-subset-files/create_selected_subset_files.sh
```

See `create_selected_subset_files.sh` and `02-subset_files.R` for more information.
The files for CI will be in `data/testing/<NEW_RELEASE>`.

### Method 3: Generate subset files for new release using a pre-existing biospecimen ID list

Use this method if you do not have a local copy of older subset files, but overall the files and samples are unchanged between releases.

To accomplish this, you will need to specify the new release to `02-subset_files.R`.
Here's an example (assuming the current directory):

```
Rscript --vanilla 02-subset_files.R \
  --biospecimen_file biospecimen_ids_for_subset.RDS \
  --output_directory ../../data/testing/<new_release> \
  --new_release <new_release>
```

Note you will still need to copy the clinical, independent sample, and BED files from the release into the testing folder, which can be done with:

```sh
FULL_DIRECTORY=data/<new_release>
SUBSET_DIRECTORY=data/testing/<new_release>

# histologies file
cp $FULL_DIRECTORY/pbta-histologies.tsv $SUBSET_DIRECTORY

# independent specimen files
cp $FULL_DIRECTORY/independent-specimens*.tsv $SUBSET_DIRECTORY

# all bed files
cp $FULL_DIRECTORY/*.bed $SUBSET_DIRECTORY

# if the md5sum.txt file already exists, get rid of it
rm -f $SUBSET_DIRECTORY/md5sum.txt
# create a new md5sum.txt file
cd $SUBSET_DIRECTORY
md5sum * > md5sum.txt

# Changelog does not get tracked
cd ../../../analyses/create-subset-files
cp $FULL_DIRECTORY/CHANGELOG.md $SUBSET_DIRECTORY
```

The files for CI will be in `data/testing/<new_release>`.

### You're ready to upload the files for CI to S3.