## Nhat Duong
## November, 22 2019

######### ASSUMPTIONS ########
# ../../scratch is available to store intermediate files
##############################

# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import subprocess
import sys
import os

# Related third party
import numpy as np
import pandas as pd

## Define the extensions
MANTA_EXT = '.manta'
CNVKIT_EXT = '.cnvkit'
FREEC_EXT = '.freec'

# Define the column headers for IDs
MANTA_ID_HEADER = 'Kids.First.Biospecimen.ID.Tumor'
CNVKIT_ID_HEADER = 'ID'
FREEC_ID_HEADER = 'Kids_First_Biospecimen_ID'

parser = argparse.ArgumentParser(description="""This script splits CNV files
                                                into one per sample. It also
                                                prints a snakemake config file
                                                to the specified filename.""")
parser.add_argument('--manta', required=True,
                    help='path to the manta file')
parser.add_argument('--cnvkit', required=True,
                    help='path to the cnvkit file')
parser.add_argument('--freec', required=True,
                    help='path to the freec file')
parser.add_argument('--snake', required=True,
                    help='path for snakemake config file')
parser.add_argument('--scratch', required=True,
                    help='directory for scratch files')
parser.add_argument('--uncalled', required=True,
                    help='path for the table of sample-caller outputs removed and not called for too many CNVs')
parser.add_argument('--maxcnvs', default=2500,
                    help='samples with more than 2500 cnvs are set to blank')
parser.add_argument('--cnvsize', default=3000,
                    help='cnv cutoff size in base pairs')
parser.add_argument('--freecp', default=0.01,
                    help='p-value cutoff for freec')

args = parser.parse_args()


## Pandas load/read files in
merged_manta = pd.read_csv(args.manta, delimiter='\t', dtype=str   )
merged_cnvkit = pd.read_csv(args.cnvkit, delimiter='\t', dtype=str )
merged_freec = pd.read_csv(args.freec, delimiter='\t', dtype=str )



## Extract the samples for each files to merge them all together. This takes into account uneven
## numbers of samples per file
manta_samples = set(merged_manta[MANTA_ID_HEADER])
cnvkit_samples = set(merged_cnvkit[CNVKIT_ID_HEADER])
freec_samples = set(merged_freec[FREEC_ID_HEADER])


## Merged and take the unique samples. Any method without a certain sample will get an empty file
## for of that sample.
all_samples = manta_samples | cnvkit_samples | freec_samples # set union

## Define and create assumed directories
scratch_d = args.scratch
manta_d = os.path.join(scratch_d, 'manta_manta')
cnvkit_d = os.path.join(scratch_d, 'cnvkit_cnvkit')
freec_d = os.path.join(scratch_d, 'freec_freec')
if not os.path.exists(manta_d):
    os.makedirs(manta_d)
if not os.path.exists(cnvkit_d):
    os.makedirs(cnvkit_d)
if not os.path.exists(freec_d):
    os.makedirs(freec_d)



bad_calls = []

## Loop through each sample, search for that sample in each of the three dataframes,
## and create a file of the sample in each directory
for sample in all_samples:

    ## Pull out the CNVs with that sample name
    manta_export = merged_manta.loc[merged_manta[MANTA_ID_HEADER] == sample]

    ## Write cnvs to file if less than maxcnvs / otherwise empty file and add to bad_calls list
    with open(os.path.join(manta_d, sample + MANTA_EXT), 'w') as file_out:
        if manta_export.shape[0] <= args.maxcnvs and manta_export.shape[0] > 0:
            manta_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tmanta\n")

    cnvkit_export = merged_cnvkit.loc[merged_cnvkit[CNVKIT_ID_HEADER] == sample]
    with open(os.path.join(cnvkit_d, sample + CNVKIT_EXT), 'w') as file_out:
        if cnvkit_export.shape[0] <= args.maxcnvs and cnvkit_export.shape[0] > 0:
            cnvkit_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tcnvkit\n")

    freec_export = merged_freec.loc[merged_freec[FREEC_ID_HEADER] == sample]
    with open(os.path.join(freec_d, sample + FREEC_EXT), 'w') as file_out:
        if freec_export.shape[0] <= args.maxcnvs and freec_export.shape[0] > 0:
            freec_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tfreec\n")


## Make the Snakemake config file. Write all of the sample names into the config file
with open(args.snake, 'w') as file:
    file.write('samples:' + '\n')
    for sample in all_samples:
        file.write('  ' + str(sample) + ':' + '\n')

    ## Define the extension for the config file
    file.write('manta_ext: ' + MANTA_EXT + '\n')
    file.write('cnvkit_ext: ' + CNVKIT_EXT + '\n')
    file.write('freec_ext: ' + FREEC_EXT + '\n')

    ## Define location for python scripts and scratch
    file.write('scripts: ' + os.path.dirname(os.path.realpath(__file__)) + '\n')
    file.write('scratch: ' + scratch_d + '\n')

    ## Define the size cutoff and freec's pval cut off.
    file.write('size_cutoff: ' + str(args.cnvsize) + '\n')
    file.write('freec_pval: ' + str(args.freecp) + '\n')

## Write out the bad calls file
with open(args.uncalled, 'w') as file:
    file.write("sample\tcaller\n")
    file.writelines(bad_calls)
