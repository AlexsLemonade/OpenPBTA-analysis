# TCGA brain tumor WXS capture kit investigation
This is an investigation of the TMB discrepancy between PBTA and TCGA data. 

## Background
To perform the analysis of [Tumor Mutation Burden Compare to TCGA](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tmb-compare), 160 tumor/normal paired TCGA WXS data had been selected and processed. The input BAM file manifest is on the `pbta-tcga-manifest.tsv` and is available in both V14 and V15 release. The MC3 WXS BED file [`gencode.v19.basic.exome.bed`](https://gdc.cancer.gov/about-data/publications/mc3-2017) was used as the capture kit region for both mutation calling and TMB calculation. It turns out the adult-based TCGA has overall lower TMB compared to pediatric-based PBTA data, and on top of that the Lancet somatic calling outputs much more low VAFs. To better figure out those issues, we started this investigation.

We started by looking only the "codon impact mutation" where we only selected the `Mutation|Splice|Silent` mutations.
- plot: https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/548#issuecomment-590813383
- pseudocode: https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/548#issuecomment-590873867

It does look like PBTA MAF has lower mutation counts than TCGA for "codon impact mutations". However the overall PBTA has more mutation counts than TCGA, that means, PBTA should have a lot more somatic mutations located in Intron and UTR regions than TCGA. We randomly picked one intron mutation from sample `BS_VVD8TJNE`, it is located at `chr1:162317128` which is the intron of gene NOS1AP, and found out it was covered by 57 reads for tumor and 21 reads for normal. Then we checked the TCGA BAM for the same loci, and found there're not reads covering that region at all. Hence we decided to investigate more on the actual capture region for TCGA data.

Genome Browser screen shot around chr1:162317128 for two TCGA samples.
![](./plots/screen-shot-genome-browser-tcga-gencode.png)

## TCGA WXS capture kit query and preparation
### 01. get-tcga-capture_kit-info
Since we got the TCGA data manifest from the [GDC portal](https://portal.gdc.cancer.gov/) by querying all the brain tumor data, and GDC actually has the capture kit name and its downloadable URL for each BAM files as part of the file metadata, so we created the python script to retrieve the actual BED files by hitting the [GDC file API endpoint](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/).

- script: [`scripts/get-tcga-capture_kit-info.py`](scripts/get-tcga-capture_kit-info.py)
- output: `results/tcga-capture_kit-info.csv`

### 02. prepare-tcga-capture_kit
By checking the `results/tcga-capture_kit-info.csv`, it turns out there are BAMs with `|` in the returned capture kit name and url which are those with more than one capture kit which neither the GDC nor its origin data center could retrieve/figure out what the actual capture kit had been applied. We plan to generate an intersected BED for those samples and used that for our analysis. We created scripts to download all unique BED files and added prefix `chr` and used [CrossMap tool](http://crossmap.sourceforge.net/) to convert all hg19 coordinates to Gh38 and saved them in  `results` folder with a `.Gh38.bed` extension

- script: 
  - [`scripts/prepare-tcga-capture_kit.sh`](scripts/prepare-tcga-capture_kit.sh)
  - [`scripts/CrossMap.sh`](scripts/CrossMap.sh)
- output: `results/*.Gh38.bed`

## Check the intersection region for the existing TCGA and PBTA MAF
The plan is to re-run the TCGA data with the new `*.Gh38.bed*` and re-do the TMB comparison using each sample's actual calling region. But to double check the new BED files, we wanted to intersect all the regions and use the overlapping region to check all the mutations counts before re-run everything.

### 03. intersect-bed-maf
We created script to prepare the input dataframe for the boxplot script. This script intersects all the BED and then use that to intersect with released PBTA and TCGA MAF and then counted all the mutation number within that intersection region, and then mapped that counts to the project(TCGA/PBTA) and tumor type for each sample.
- script: [`scripts/intersect-bed-maf.sh`](scripts/intersect-bed-maf.sh)
- output: `scratch/somatic-count_with-histologies.tsv`

### 04. mutation-counts-boxplot
- script: [`scripts/boxplot.R`](scripts/boxplot.R)
- output: `plots/boxplot-*.png`

![](plots/boxplot-all.png)

## Coverage-Comparison
In order to show in the form of coverage that the downloaded BED files were  correct, we calculated the number of bases that are covered  at least at 20x in both old TCGA BED file(`gencode.v19.basic.exome`) and the new downloaded BED files (`results/*Gh38.bed`). The following command was used to run sambamba coverage per base with a coverage threshold of 20x -
`sambamba_v0.5.9  depth base -c 20 -L *Gh38.bed *.bam > *.sambamba_basecoverage.txt`
The coverage table in `results/TCGA_oldandnew_coverage_comparisons.txt` gives  a  tabular summary of  number of bases and percentages within corresponding BED files. Boxplot comparison for the same data is available here(`plots/TCGA_oldandnew_coverage_plots.png`)

![](plots/TCGA_oldandnew_coverage_plots.png)

## (Planned) rerun-tcga-with-new-bed
