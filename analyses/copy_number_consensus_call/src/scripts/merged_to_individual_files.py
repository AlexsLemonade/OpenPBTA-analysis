## Nhat Duong
## November, 22 2019

import numpy as np
import pandas as pd
import subprocess
import sys
import os

######### ASSUMPTIONS ########
# 1) Files are feed into stdout as manta, cnvkit, then freec
# 2) That the ../../scratch directory are always available to store intermediate files
# 3) Single sample with more than 2500 CNVs are not taken into account because that's a lot of CNVs. Empty file will be created
# 4) CNV size cut-off is 3000 bp
# 5) pval to filter out for freec is 0.01


## Get the list of file names from stdin
list_of_files = []
for a in sys.argv[1:]:
    list_of_files += [a]

## Check if there are 3 files, if not, abort mission
if len(list_of_files) != 3:
    sys.stderr.write('Need 3 input files for the code to work on. Please check and retry\n')
    sys.exit()


## Define the extension based on the order of manta, cnvkit, and freec
manta_ext = '.manta'
cnvkit_ext = '.cnvkit'
freec_ext = '.freec'


## Assign the files from stdin
manta_gz_file = list_of_files[0]
cnvkit_gz_file = list_of_files[1]
freec_gz_file = list_of_files[2]

## Pandas load/read files in
merged_manta = pd.read_csv(manta_gz_file,delimiter='\t')
merged_cnvkit = pd.read_csv(cnvkit_gz_file,delimiter='\t')
merged_freec = pd.read_csv(freec_gz_file,delimiter='\t')


## Extract the samples for each files to merge them all together. This takes into account uneven
## numbers of samples per file
manta_samples = np.unique(merged_manta['Kids.First.Biospecimen.ID.Tumor'])
cnvkit_samples = np.unique(merged_cnvkit['ID'])
freec_samples = np.unique(merged_freec['Kids_First_Biospecimen_ID'])


## Merged and take the unique samples. Any method without a certain sample will get an empty file
## for of that sample.
all_samples = np.unique(list(manta_samples) + list(cnvkit_samples)  + list(freec_samples))

## Define and create assumed directories
manta_d = '../../scratch/manta_manta'
cnvkit_d = '../../scratch/cnvkit_cnvkit'
freec_d = '../../scratch/freec_freec'
if not os.path.exists(manta_d):
    os.makedirs(manta_d)
if not os.path.exists(cnvkit_d):
    os.makedirs(cnvkit_d)
if not os.path.exists(freec_d):
    os.makedirs(freec_d)


## Loop through each sample, search for that sample in each of the three dataframes,
## and create a file of the sample in each directory
for sample in all_samples:

    ## Pull out the CNVs with that sample name
    manta_export = merged_manta.loc[merged_manta['Kids.First.Biospecimen.ID.Tumor'] == sample]

    ## Open up a file and write to it if the number of CNVs are less than 2500
    with open((manta_d + '/' + sample + manta_ext), 'w') as file_out:
        if manta_export.shape[0] <= 2500:
            manta_export.to_csv(file_out, sep='\t', index=False)
        else:
            pass

    cnvkit_export = merged_cnvkit.loc[merged_cnvkit['ID'] == sample]
    with open((cnvkit_d + '/' + sample + cnvkit_ext), 'w') as file_out:
        if cnvkit_export.shape[0] <= 2500:
            cnvkit_export.to_csv(file_out, sep='\t', index=False)
        else:
            pass

    freec_export = merged_freec.loc[merged_freec['Kids_First_Biospecimen_ID'] == sample]
    with open((freec_d + '/' + sample + freec_ext), 'w') as file_out:
        if freec_export.shape[0] <= 2500:
            freec_export.to_csv(file_out, sep='\t', index=False)
        else:
            pass


## Make the Snakemake config file. Write all of the sample names into the config file
with open('../../scratch/config_snakemake.yaml', 'w') as file:
    file.write('samples:' + '\n')
    for sample in all_samples:
        file.write('  ' + str(sample) + ':' + '\n')

## Define the extension for the config file
    file.write('manta_ext: ' + manta_ext + '\n')
    file.write('cnvkit_ext: ' + cnvkit_ext + '\n')
    file.write('freec_ext: ' + freec_ext + '\n')

## Define location for python scripts
    file.write('scripts: src/scripts/' + '\n')

## Define the size cutoff and freec's pval cut off.
    file.write('size_cutoff: 3000' + '\n')
    file.write('freec_pval: 0.01' + '\n')
