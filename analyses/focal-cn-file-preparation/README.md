## Focal copy number file preparation

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

The copy number data from OpenPBTA are provided as ranges or segments.
The purpose of this module is to map from those ranges to gene identifiers for consumption by downstream analyses (e.g., OncoPrint plotting).

### Running this analysis

To run this analysis, use the following (from the root directory of the repository):

```
bash analyses/focal-cn-file-preparation/run-prepare-cn.sh
```

### Scripts and notebooks

* `01-add-ploidy-cnvkit.Rmd` - The two CNV callers, CNVkit and ControlFreeC, do not handle ploidy in the same way ([A Note on Ploidy](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/de661fbe740717472fcf01c7d9b74fe1b946aece/doc/data-formats.md#a-note-on-ploidy) in the Data Formats documentation). 
  This notebook adds the ploidy inferred via ControlFreeC to the CNVkit data and adds a status column that defines gain and loss broadly.
  (Note that [the logic around sex chromosomes in males when ploidy = 3 leaves something to be desired](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/259#discussion_r345354403)).
* `02-add-ploidy-consensus.Rmd` - This is very similar to the CNVkit file prep.
However, there are instances in the consensus SEG file where `copy.num` is `NA` which are removed.
See the notebook for more information.
* `03-prepare-cn-file.R` - This script performs the ranges to annotation mapping using the GENCODE v27 GTF included via the data download step; it takes the ControlFreeC file or a SEG (e.g., CNVkit, consensus SEG) file prepared with `01-add-ploidy-cnvkit.Rmd` and  `02-add-ploidy-cnvkit.Rmd`as input.
  The mapping is limited to _exons_. 
  Mapping to cytobands is performed with the [`org.Hs.eg.db`](https://doi.org/doi:10.18129/B9.bioc.org.Hs.eg.db) package.
  A table with the following columns is returned:
  
  | biospecimen_id | status | copy_number | ploidy | ensembl | gene_symbol | cytoband |
  |----------------|--------|-------------|--------|---------|-------------|---------|
  
* `rna-expression-validation.R` - This script examines RNA-seq expression levels (RSEM FPKM) of genes that are called as deletions. 
It is not currently run via the shell script.

### Output files for downstream consumption 
  
```
results
├── cnvkit_annotated_cn_autosomes.tsv.gz
├── cnvkit_annotated_cn_x_and_y.tsv.gz
├── consensus_seg_annotated_cn_autosomes.tsv.gz
├── consensus_seg_annotated_cn_x_and_y.tsv.gz
├── controlfreec_annotated_cn_autosomes.tsv.gz
└── controlfreec_annotated_cn_x_and_y.tsv.gz
```
