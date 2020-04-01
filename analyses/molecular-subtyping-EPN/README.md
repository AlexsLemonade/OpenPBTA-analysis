# Molecular Subtyping of Ependymoma 

<b>Module Authors:</b> Teja Koganti(<a href="https://github.com/tkoganti">@tkoganti</a>) and Josh Shapiro(<a href="https://github.com/jashapiro">@jashapiro</a>)

In this analysis we subtype ependymoma samples based on fusions, CNV, NFKB_pathway_GSEAscore, breaks_density-chromosomal_instability_CNV, breaks_density-chromosomal_instability_SV, GISTIC_focal_CN_CDKN2A and gene expression data

## Usage
`bash run-molecular-subtyping-EPN.sh`

This above  script is designed to change to this directory to run, so it should run from any location.

## Folder content

1. <b>`00-subset-for-EPN.R`</b> is a script that takes subsets only  ependymoma samples expression data for CI. The script uses `pbta-histologies.tsv` file to filter for ependymoma samples and     `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` file for expression  data

2. <b>`01-make_notebook_RNAandDNA.py`</b> is a script that filters WGS  and RNA-seq ependymoma samples from `pbta-histologies.tsv` file and adds `disease_group` column to the output file based on the primary_site. The values for `disease_group` are supratentorial/infratentorial/undetermined

3. <b>`02_ependymoma_generate_all_data.py`</b>  is a script that takes in expression, GISTIC, fusion, breakpoint, GISTIC, GSVA files to add values from these tables as new columns to the input notebook. Output from `01-make_notebook_RNAandDNA.py` script is used as input notebook. The output notebook from this is saved to `results/EPN_all_data.tsv`

4. <b> `03-subgrouping_samples.py`  </b>  is a script that takes the table `results/EPN_all_data.tsv`  as input and adds a column that groups the samples into one of these groups - ST-EPN-RELA, ST-EPN-YAP1, PF-EPN-A, and PF-EPN-B
