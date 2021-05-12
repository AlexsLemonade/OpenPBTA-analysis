# Copy Number Consensus Call

## Overview

The PBTA data set contains CNVs called from different callers, ie. Manta, CNVkit, and Freec. 
The goal is to use all of these callers to reduce false positives and come up with a final consensus list of CNVs.
This analysis uses information from the following files generated from the 3 callers

* `pbta-cnv-cnvkit.seg.gz`
* `pbta-cnv-controlfreec.tsv.gz`
* `pbta-sv-manta.tsv.gz`

The analysis produces the following output files

* `results/cnv_consensus.tsv`:  A tab separated file of consensus copy number variants, including the original calls used for each consensus call.
  *Note: only samples that pass QC are included, see below*
* `results/uncalled_samples.tsv`: A table of sample-caller pairs excluded from the consensus due to QC failure (usually too many CNV calls, see [Methods: Consensus CNV creation](#consensus-cnv-creation)) If a sample fails QC from two or more callers, it will not appear in downstream files.
* `results/pbta-cnv-consensus.seg.gz`: A `.seg` formatted file for downstream processing. *Only includes samples that pass QC*

* `ref/cnv_excluded_regions.bed`: A `.bed` file of error-prone regions that were filtered from copy number calls
* `ref/cnv_callable.bed`: A `.bed` file of regions considered "callable" by the analysis pipeline

## Running the pipeline

To run the entire pipeline, make sure to have the latest release of the three input files mentioned in the Overview section.
Go to OpenPBTA-analysis/analyses/copy_number_consensus_call and run `bash run_consensus_call.sh`

## Methods

### Assayed Regions

Regions of the genome with a high potential for error are first defined by merging the set telomeric, centromeric and heterochromatic regions with regions around immunoglobulins and segmentmental duplications.
The input files for this step are described in `scripts/prepare_blacklist_files.sh` and include:

* `ref/centromeres.bed`
* `ref/heterochromatin.bed`
* `ref/immunoglobulin_regions.bed`
* `ref/segmental_dups.bed`
* `ref/telomeres.bed`

The final set of merged excluded regions are placed in the file `ref/cnv_excluded_regions.bed`

In addition, a file of the genomic regions that we deem "callable" is created at `ref/cnv_callable.bed` as the complement of the excluded regions, after removing exclusions smaller than 200kb.

### Consensus CNV creation

The per-sample pipeline revolves around the use of Snakemake to run analysis for each patient sample. The overview of the steps are as followed:

1) Parse through the 3 input files and put CNVs of the **same caller and sample** in the same files.
2) Remove any sample/caller combination files with **more than 2500** CNVs called.
   We belive these to be noisy/poor quality samples (this came from what GISTIC uses as a cutoff for noisy samples).
3) Create a `config_snakemake.yaml` that contains all of the samples names to run the Snakemake pipeline
4) Run the Snakemake pipeline to perform analysis **per sample**. 
5) Filter for any CNVs that are over a certain **SIZE_CUTOFF** (default 3000 bp)
6) Filter for any **significant** CNVs called by Freec (default pval = 0.01)
7) Filter out any CNVs that overlap 50% or more with **Immunoglobulin, telomeric, centromeric, seg_dup regions** as found in the file `ref/cnv_excluded.bed`
8) Merge any CNVs of the same sample and call method if they **overlap or within 10,000 bp** (We consider CNV calls within 10,000 bp the same CNV)
9) Reformat the columns of the files (So the info are easier to read)
10) **Call consensus** by comparing CNVs from 2 call methods at a time. 

Since there are 3 callers, there were 3 comparisons: `manta-cnvkit`, `manta-freec`, and `cnvkit-freec`. If a CNV from 1 caller **overlaps 50% or more** with at least 1 CNV from another caller, the common region of the overlapping CNV would be the new CONSENSUS CNV.

11) **Sort and merge** the CNVs from the comparison pairs ,`manta-cnvkit` `manta-freec` `cnvkit-freec`, together into 1 file
12) Resolve overlapping segments where duplications are embedded within larger deletion segments, or deletions within duplications.
13) After every samples' consensus CNVs were called, **combine all merged files** from step 10 and output to `results/cnv_consensus.tsv`
14) The `results/cnv_consensus.tsv` is translated into a `results/pbta-cnv-consensus.seg` file in the same format as `pbta-cnv-cnvkit.seg.gz`, including all samples where at least two callers passed quality filtering.
When a consensus segment is derived from multiple source segments, we take the mean of the CNVkit `seg.mean` values from the source segments, weighted by segment length.
If no CNVkit variant was included, the value for this column is NA.
The `copy.num` column is the weighted median of Control-FREEC segment values where they exist, or CNVkit values in the absence of Control-FREEC data.
Because some software (notably GISTIC) requires all samples to have the same regions called, the copy number variants from `cnv_consensus.tsv` are supplementented with "neutral" segments where no call was made.
These include all non-variant regions present in `ref/cnv_callable.bed`
The neutral regions are assigned NA.

## Example Output File

```
chrom	start	end	manta_CNVs	cnvkit_CNVs	freec_CNVs	CNV_type	Biospecimen
chr11	771036	866778	NULL	770516:866778:3:0.214821	771036:871536:3:NA	DUP	BS_007JTNB8
chr13	99966948	99991872	NULL	99954829:99994557:3:0.263969	99966948:99991872:3:NA	DUP	BS_007JTNB8
chr14	103515996	103563240	NULL	103515996:103563363:3:0.237098	103511784:103541532:3:NA,103543140:103563240:3:NA	DUP	BS_007JTNB8
```
* The 1st line of the file is the header which contains the column names. There are 9 columns in total
* The 2nd line is the first CNV of the file.
* Column 1 is the **consensus** CNV chromosome
* Column 2 is the **consensus** CNV start location
* Column 3 is the **consensus** CNV end location
* Columns 4, 5, and 6 contain the calls of Manta, CNVkit, and Freec that make up the **consensus** CNV described in columns 1, 2, and 3. 
* ie. If there is info in column 4, that means one or more CNVs called from Manta made up the current **consensus** CNV described in columns 1, 2, and 3. 
* Columns 4, 5, and 6 have the following format: `START:END:COPY_NUMBER,START:END:COPY_NUMBER:SEG.MEAN`
  * Note that if there is more than one original CNV call corresponding to a given consensus CNV from a given caller, the information for each of the CNV calls will be comma separated.
  * In the example output above column 6 of line 4 contains `103511784:103541532:3:NA,103543140:103563240:3:NA` which means 2 CNVs called by FreeC helped to make up the **consensus** CNV on line 4. 
One has the start and end coordinates of `103511784:103541532` **on the same chromosome**, has a copy number of `3`, and a seg.mean of `NA`. The other has the coordinates `103543140:103563240`, has a copy number of `3`, and a seg.mean of `NA`
* Column 7 is the CNVtype. This will be one of DUP or DEL, corresponding to duplications or deletions, respectively. Note that this does not describe the number of copies, only the direction of the copy number change.
* Column 8 is the Biospecimen Sample name
