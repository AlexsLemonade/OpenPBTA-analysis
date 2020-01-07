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



## Running the pipeline

To run the entire pipeline, make sure to have the latest release of the three input files mentioned in the Overview section.
Go to OpenPBTA-analysis/analyses/copy_number_consensus_call and run `bash run_consensus_call.sh`



