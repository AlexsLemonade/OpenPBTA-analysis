## Focal copy number file preparation

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)), Joshua Shapiro ([@jashapiro](https://github.com/jashapiro)), and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

The copy number data from OpenPBTA are provided as ranges or segments.
The purpose of this module is to map from those ranges to gene identifiers for consumption by downstream analyses (e.g., OncoPrint plotting).

### Running this analysis
*This analysis requires at least ~24 GB of RAM to run to completion*

To run this analysis _only on consensus SEG file_, 

use OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder and relative path to `copy_number_consensus_call/results/pbta-cnv-consensus.seg.gz` while running molecular-subtyping modules for release.

```
OPENPBTA_BASE_SUBTYPING=1 bash analyses/focal-cn-file-preparation/run-prepare-cn.sh
```

Or by default runs analyses using pbta-histologies.tsv and downloaded files from data release:

```
bash analyses/focal-cn-file-preparation/run-prepare-cn.sh
```

**Note**: The `run-bedtools.snakemake` script is implemented in `run-prepare-cn.sh` to run the bedtools coverage steps between the UCSC cytoband file and the samples in the copy number files produced in `02-add-ploidy-consensus.Rmd`.
This script currently takes a while to run, and therefore slows down the processing speed of the main shell script `run-prepare-cn.sh`.

Running the following from the root directory of the repository will run the steps for the original copy number call files (CNVkit and ControlFreeC) in addition to the consensus SEG file:

```
RUN_ORIGINAL=1 bash analyses/focal-cn-file-preparation/run-prepare-cn.sh
```

### Scripts and notebooks

* `01-add-ploidy-cnvkit.Rmd` - The two CNV callers, CNVkit and ControlFreeC, do not handle ploidy in the same way ([A Note on Ploidy](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/de661fbe740717472fcf01c7d9b74fe1b946aece/doc/data-formats.md#a-note-on-ploidy) in the Data Formats documentation).
  This notebook adds the ploidy inferred via ControlFreeC to the CNVkit data and adds a status column that defines gain and loss broadly.
  Specifically, segments with copy number fewer than ploidy are losses, segments with copy number greater than ploidy are marked as a gain, and segments where copy number is equal to ploidy are marked as neutral.
  (Note that [the logic around sex chromosomes in males when ploidy = 3 leaves something to be desired](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/259#discussion_r345354403)).

* `02-add-ploidy-consensus.Rmd` - This is very similar to the CNVkit file prep (`01-add-ploidy-cnvkit.Rmd`).
However, there are instances in the consensus SEG file where `copy.num` is `NA` which are removed.
See the notebook for more information. This notebook also prepares lists of copy number bed files by sample for use in the implementation of bedtools coverage in `run-bedtools.snakemake`.

* `03-add-cytoband-status-consensus.Rmd` - This notebook reads in the bedtools coverage output files and defines the dominant copy number status for each chromosome arm. The output of this notebook is a table with the following columns:

  | `Kids_First_Biospecimen_ID` | chr | cytoband | dominant_status | band_length | callable_fraction | gain_fraction | loss_fraction | chromosome_arm |
  |----------------|--------|-------------|--------|---------|----------|-------------|---------|---------------|

* `04-prepare-cn-file.R` - This script performs the ranges to annotation mapping using the GENCODE v27 GTF included via the data download step; it takes the ControlFreeC file or a SEG (e.g., CNVkit, consensus SEG) file prepared with `01-add-ploidy-cnvkit.Rmd` and  `02-add-ploidy-cnvkit.Rmd` as input.
  **The mapping is limited to _exons_.**
  Mapping to cytobands is performed with the [`org.Hs.eg.db`](https://doi.org/doi:10.18129/B9.bioc.org.Hs.eg.db) package.
  A table with the following columns is returned:

  | biospecimen_id | status | copy_number | ploidy | ensembl | gene_symbol | cytoband |
  |----------------|--------|-------------|--------|---------|-------------|---------|
  Any segment that is copy neutral is filtered out of this table. In addition, [any segments with copy number > (2 * ploidy) are marked as amplifications](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/e2058dd43d9b1dd41b609e0c3429c72f79ff3be6/analyses/focal-cn-file-preparation/03-prepare-cn-file.R#L275) in the `status` column. 
Any segments with copy number == 0  are marked as deep deletions in the `status` column.

* `05-define-most-focal-cn-units.Rmd` - This notebook defines the _most focal_ recurrent copy number units by removing focal changes that are within entire chromosome arm losses and gains.
_Most focal_ here meaning if a chromosome arm is not clearly defined as a gain or loss (and is callable) we look to define the cytoband level status.
Similarly, if a cytoband is not clearly defined as a gain or loss (and is callable) we then look to define the gene level status.
To make these calls, the following decisions around cutoffs were made:
	- The percentage of calls a particular status needs to be above to be called the majority status is currently 	90%.
	This decision was made to ensure that the status is not only the majority status but is also 	significantly called more than the other status values in the region.
	- The percentage of a region (arm, cytoband, or gene) that should be callable is more than 50%.
	This decision was made because it seems reasonable to expect a region to be more than 50% callable for a 	dominant status call to be made.

* `06-find-recurrent-calls.Rmd` - This notebook determines the recurrent focal copy number dominant status calls by region using the output of `05-define-most-focal-cn-units.Rmd`.
Recurrence here has been arbitrarily defined based on the plotting of the distribution of status calls and a [similar decision](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/66bb67a7bf29aad4510a0913a2dbc88da0013be8/analyses/fusion_filtering/06-recurrent-fusions-per-histology.R#L152) made in `analyses/fusion_filtering/06-recurrent-fusions-per-histology.R` to make the cutoff for recurrence to be greater than a count of 3 samples that have the same CN status call in the same region.
This notebook returns a `TSV` file with the recurrent copy number status calls, regions and biospecimen IDs.

* `rna-expression-validation.R` - This script examines RNA-seq expression levels (RSEM FPKM) of genes that are called as deletions.
It produces loss/neutral and zero/neutral correlation plots, as well as stacked barplots displaying the distribution of ranges in expression across each of the calls (loss, neutral, zero).
_Note: The shell script's default behavior is to produce these plots using the annotated consensus SEG autosome and sex chromsome files found in this module's `results` directory and listed below._


### Output files for downstream consumption

**Note:** The output files from `03-prepare-cn-file.R` have neutral calls filtered out to reduce file size.

```
results
├── cnvkit_annotated_cn_autosomes.tsv.gz
├── cnvkit_annotated_cn_x_and_y.tsv.gz
├── consensus_seg_annotated_cn_autosomes.tsv.gz
├── consensus_seg_annotated_cn_x_and_y.tsv.gz
├── consensus_seg_most_focal_cn_status.tsv.gz
├── consensus_seg_recurrent_focal_cn_units.tsv
├── consensus_seg_with_ucsc_cytoband_status.tsv.gz
├── controlfreec_annotated_cn_autosomes.tsv.gz
└── controlfreec_annotated_cn_x_and_y.tsv.gz
```

### Folder Structure

```
focal-cn-file-preparation
├── 01-add-ploidy-cnvkit.Rmd
├── 01-add-ploidy-cnvkit.nb.html
├── 02-add-ploidy-consensus.Rmd
├── 02-add-ploidy-consensus.nb.html
├── 03-add-cytoband-status-consensus.Rmd
├── 03-add-cytoband-status-consensus.nb.html
├── 04-prepare-cn-file.R
├── 05-define-most-focal-cn-units.Rmd
├── 05-define-most-focal-cn-units.nb.html
├── 06-find-recurrent-calls.Rmd
├── 06-find-recurrent-calls.nb.html
├── README.md
├── annotation_files
│   └── txdb_from_gencode.v27.gtf.db
├── display-plots.md
├── driver-lists
├── gistic-results
│   └── pbta-cnv-cnvkit-gistic
│       ├── D.cap1.5.mat
│       ├── all_data_by_genes.txt
│       ├── all_lesions.conf_90.txt
│       ├── all_thresholded.by_genes.txt
│       ├── amp_genes.conf_90.txt
│       ├── amp_qplot.pdf
│       ├── amp_qplot.png
│       ├── broad_data_by_genes.txt
│       ├── broad_gistic_plot.pdf
│       ├── broad_significance_results.txt
│       ├── broad_values_by_arm.txt
│       ├── del_genes.conf_90.txt
│       ├── del_qplot.pdf
│       ├── del_qplot.png
│       ├── focal_dat.0.98.mat
│       ├── focal_data_by_genes.txt
│       ├── freqarms_vs_ngenes.pdf
│       ├── gistic_inputs.mat
│       ├── peak_regs.mat
│       ├── perm_ads.mat
│       ├── raw_copy_number.pdf
│       ├── raw_copy_number.png
│       ├── regions_track.conf_90.bed
│       ├── sample_cutoffs.txt
│       ├── sample_seg_counts.txt
│       ├── scores.0.98.mat
│       ├── scores.gistic
│       └── wide_peak_regs.mat
├── plots
│   ├── cnvkit_annotated_cn_autosomes_polya_loss_cor_plot.png
│   ├── cnvkit_annotated_cn_autosomes_polya_stacked_plot.png
│   ├── cnvkit_annotated_cn_autosomes_polya_zero_cor_plot.png
│   ├── cnvkit_annotated_cn_autosomes_stranded_loss_cor_plot.png
│   ├── cnvkit_annotated_cn_autosomes_stranded_stacked_plot.png
│   ├── cnvkit_annotated_cn_autosomes_stranded_zero_cor_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_polya_loss_cor_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_polya_stacked_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_polya_zero_cor_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_stranded_loss_cor_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_stranded_stacked_plot.png
│   ├── cnvkit_annotated_cn_x_and_y_stranded_zero_cor_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_polya_loss_cor_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_polya_stacked_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_polya_zero_cor_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_stranded_loss_cor_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_stranded_stacked_plot.png
│   ├── consensus_seg_annotated_cn_autosomes_stranded_zero_cor_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_polya_loss_cor_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_polya_stacked_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_polya_zero_cor_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_stranded_loss_cor_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_stranded_stacked_plot.png
│   ├── consensus_seg_annotated_cn_x_and_y_stranded_zero_cor_plot.png
│   ├── controlfreec_annotated_cn_autosomes_polya_loss_cor_plot.png
│   ├── controlfreec_annotated_cn_autosomes_polya_stacked_plot.png
│   ├── controlfreec_annotated_cn_autosomes_polya_zero_cor_plot.png
│   ├── controlfreec_annotated_cn_autosomes_stranded_loss_cor_plot.png
│   ├── controlfreec_annotated_cn_autosomes_stranded_stacked_plot.png
│   ├── controlfreec_annotated_cn_autosomes_stranded_zero_cor_plot.png
│   ├── controlfreec_annotated_cn_x_and_y_polya_loss_cor_plot.png
│   ├── controlfreec_annotated_cn_x_and_y_polya_stacked_plot.png
│   ├── controlfreec_annotated_cn_x_and_y_polya_zero_cor_plot.png
│   ├── controlfreec_annotated_cn_x_and_y_stranded_loss_cor_plot.png
│   ├── controlfreec_annotated_cn_x_and_y_stranded_stacked_plot.png
│   └── controlfreec_annotated_cn_x_and_y_stranded_zero_cor_plot.png
├── results
│   ├── cnvkit_annotated_cn_autosomes.tsv.gz
│   ├── cnvkit_annotated_cn_x_and_y.tsv.gz
│   ├── consensus_seg_annotated_cn_autosomes.tsv.gz
│   ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
│   ├── consensus_seg_most_focal_cn_status.tsv.gz
│   ├── consensus_seg_recurrent_focal_cn_units.tsv
│   ├── consensus_seg_with_ucsc_cytoband_status.tsv.gz
│   ├── controlfreec_annotated_cn_autosomes.tsv.gz
│   └── controlfreec_annotated_cn_x_and_y.tsv.gz
├── rna-expression-validation.R
├── run-bedtools.snakemake
├── run-prepare-cn.sh
└── util
    └── rna-expression-functions.R
```
