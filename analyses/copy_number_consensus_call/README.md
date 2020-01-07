# Copy Number Consensus Call

## Overview

The PBTA data set contains CNVs called from different callers, ie. Manta, CNVkit, and Freec. 
The goal is to use all of these callers to reduce false positives and come up with a final consensus list of CNVs.
This analysis uses information from the following files generated from the 3 callers

`pbta-cnv-cnvkit.seg.gz`	

`pbta-cnv-controlfreec.tsv.gz`

`pbta-sv-manta.tsv.gz`


## Methods

This pipeline revolves around the use of Snakemake to run analysis for each patient sample. The overview of the steps are as followed:

1) Parse through the 3 input files and put CNVs of the **same caller and sample** in the same files.
2) Get rid of any one sample with **more than 2500** CNVs called
2) Create a `config_snakemake.yaml` that contains all of the samples names to run the Snakemake pipeline
3) Run the Snakemake pipeline to perform analysis **per sample**. 
4) Filter for any CNVs that are over a certain **SIZE_CUTOFF** (default 3000 bp)
5) Filter for any **significant** CNVs called by Freec (default pval = 0.01)
6) Filter out any CNVs that overlap 50% or more with **IGLL, telomeric, centromeric, seg_dup regions**
7) Merge any CNVs of the same sample and call method if they **overlap or within 10,000 bp** (We consider CNV calls within 10,000 bp the same CNV)
8) Reformat the columns of the files (So the info are easier to read)
9) **Call consensus** by comparing CNVs from 2 call methods at a time. Since there are 3 callers, there were 3 comparisons
`manta-cnvkit` `manta-freec` `cnvkit-freec`. If a CNV from 1 caller **overlaps 50% or more** with at least 1 CNV from another caller,
the common region of the overlapping CNV would be the new CONSENSUS CNV. 
10) **Sort and merge** the CNVs from the comparison pairs ,`manta-cnvkit` `manta-freec` `cnvkit-freec`, together into 1 file
11) After every samples' consensus CNVs were called, **combine every single files** from step 10 and output it to `results/cnv_consensus.tsv`

## Example Output File

```
chrom start end	manta_CNVs	cnvkit_CNVs	freec_CNVs	CNV_type	sample	file_names
chr11 771036  866778	NULL	770516:866778:3	771036:871536:3	DUP	BS_007JTNB8	BS_007JTNB8.cnvkit_freec.dup.bed
chr13	99966948	99991872	NULL	99954829:99994557:3	99966948:99991872:3	DUP	BS_007JTNB8	BS_007JTNB8.cnvkit_freec.dup.bed
chr14	103515996	103563240	NULL	103515996:103563363:3	103511784:103541532:3,103543140:103563240:3	DUP	BS_007JTNB8	BS_007JTNB8.cnvkit_freec.dup.bed
```

* The 1st line of the file is the header which contains the column names. There are 9 columns in total
* The 2nd line is the first CNV of the file.
* Column 1 is the **consensus** CNV chromosome
* Column 2 is the **consensus** CNV start location
* Column 3 is the **consensus** CNV end location
* Columns 4, 5, and 6 contain the calls of Manta, CNVkit, and Freec that make up the **consensus** CNV described in columns 1, 2, and 3. 
* ie. If there is info in column 4, that means one or more CNVs called from Manta made up the current **consensus** CNV described in columns 1, 2, and 3. 
* Columns 4, 5, and 6 have the following format `START:END:COPY_NUMBERS,START:END:COPY_NUMBERS`
* ie. Take a look at the code block above of the output file. Column 6 of line 4 contains `103511784:103541532:3,103543140:103563240:3` which means 2 CNVs called by FreeC helped to make up the **consensus** CNV on line 4. One has the coordinate of `103511784:103541532` **on the same chromosome** and has a copy number of 3 and another one has the coordinate of `103543140:103563240` **on the same chromosome** and has a copy number of 3. 
* Column 7 is the CNVtype. In this case, it is either DUP or DEL
* Column 8 is the Sample name
* Column 9 contains the name of of the files (`manta-cnvkit` `manta-freec` `cnvkit-freec`) that made up the **consensus** CNV. 



## Running the pipeline

To run the entire pipeline, make sure to have the latest release of the three input files mentioned in the Overview section.
Go to OpenPBTA-analysis/analyses/copy_number_consensus_call and run `bash run_consensus_call.sh`



