## Running the figure generation script

Steps to refresh all publication-ready figures.
All steps assume your current directory is the top of this repository.

#### 1. Obtain the current dataset.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#how-to-obtain-openpbta-data) to obtain the current data release.
We recommend [using the download script](https://github.com/AlexsLemonade/OpenPBTA-analysis#data-access-via-download-script) to obtain data because this will automatically create symlinks in `data/` to the latest files.

#### 2. Set up an up-to-date project Docker container.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-image) for setting up the project Docker container.
Briefly, the latest version of the project Docker image, which is updated upon commit to `master`, can be obtained and run via:
```
docker pull ccdlopenpbta/open-pbta:latest
docker run \
  -e PASSWORD=<password> \
  -p 8787:8787 \
  -v $(pwd):/home/rstudio/kitematic \
  ccdlopenpbta/open-pbta:latest
```
You may choose to use [`docker exec`](https://docs.docker.com/engine/reference/commandline/exec/) to interact with the container from there or if you'd prefer the RStudio interface, you can navigate to `localhost:8787` and enter username `rstudio` and the password you set with the `run` command above.

#### 3. Run the bash script that generates the figures (`scripts/run-figures.sh`).

This script runs **_all_** the intermediate steps needed to generate figures starting with the original data files.

```
bash figures/generate-figures.sh
```

Figures are saved to the `figures/pngs` folder and will be linked to the accompanying manuscript repository [`AlexsLemonade/OpenPBTA-manuscript`](https://github.com/AlexsLemonade/OpenPBTA-manuscript/).

## Summary for each figure

Each figure has its own script stored in the `figures/scripts`.
All are called by the main bash script `figures/run-figures.sh`.
However, we list information about the resources, intermediate steps, and [PBTA data files](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/data-formats.md#pbta-data-files) required for generating each figure below for convenience.

| Figure | Individual script | Notes on requirements | Linked analysis modules | PBTA data files consumed |
|--------|--------|------------------|-------------------------|-----------------------------|
| Figure 1 | No individual script | No high RAM requirements | [`sample-distribution-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/sample-distribution-analysis) |  `pbta-histologies.tsv` |
| Figure 2 | [`scripts/fig2-mutational-landscape.R`](./scripts/fig2-mutational-landscape.R) | 256GB of RAM are needed due to the run_caller_consensus_analysis-pbta.sh handling of large MAF files|[`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) <br> [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures) |  `pbta-snv-lancet.vep.maf.gz` <br> `pbta-snv-mutect2.vep.maf.gz` <br> `pbta-snv-strelka2.vep.maf.gz` <br> `pbta-snv-vardict.vep.maf.gz` <br> `tcga-snv-lancet.vep.maf.gz` <br> `tcga-snv-mutect2.vep.maf.gz` <br> `tcga-snv-strelka2.vep.maf.gz` |
| CN status heatmap | [`analyses/copy_number_consensus_call/run_consensus_call.sh`](./analyses/copy_number_consensus_call/run_consensus_call.sh) and [`analyses/cnv-chrom-plot/cn_status_heatmap.Rmd`](./analyses/cnv-chrom-plot/cn_status_heatmap.Rmd) | No high RAM requirements | [`cnv-chrom-plot`](./analyses/cnv-chrom-plot) |  `pbta-cnv-controlfreec.tsv.gz` <br> `pbta-sv-manta.tsv.gz` <br> `pbta-cnv-cnvkit.seg.gz` |
| Figure 3 | No individual script <br> ([`analyses/focal-cn-file-preparation/run-prepare-cn.sh`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/focal-cn-file-preparation/run-prepare-cn.sh) and [`analyses/oncoprint-landscape/run-oncoprint.sh`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/oncoprint-landscape/run-oncoprint.sh) scripts are used)              | 24GB of RAM are needed due to the `run-prepare-cn.sh` handling of large copy number files | [`focal-cn-file-preparation`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/focal-cn-file-preparation) <br> [`oncoprint-landscape`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/oncoprint-landscape) | `pbta-histologies.tsv` <br> `pbta-snv-consensus-mutation.maf.tsv.gz` <br> `pbta-fusion-putative-oncogenic.tsv` <br> `consensus_seg_annotated_cn_autosomes.tsv.gz` <br> `independent-specimens.wgs.primary-plus.tsv` |
| Transcriptomic overview | [scripts/transcriptomic-overview.R](./scripts/transcriptomic-overiew.R) | Due to the GSVA steps, we recommend ~32 GB of RAM for generating this figure | [`transcriptomic-dimension-reduction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/transcriptomic-dimension-reduction) <br> [`collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) <br> [`gene-set-enrichment-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/gene-set-enrichment-analysis) <br> [`immune-deconv`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv) | `pbta-histologies.tsv` <br> `pbta-gene-expression-rsem-fpkm.stranded.rds` |
| Mutation co-occurrence | No individual script <br> ([`analyses/interaction-plots/01-create-interaction-plots.sh`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/interaction-plots/01-create-interaction-plots.sh) is used) |  No high RAM requirements | [`interaction-plots`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/interaction-plots) |`independent-specimens.wgs.primary-plus.tsv` <br> `pbta-snv-consensus-mutation.maf.tsv.gz`  |
| Telomerase activities | [scripts/TelomeraseActivitites.R](./scripts/TelomeraseActivitites.R)|  No high RAM requirements |[`telomerase-activity-prediction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction/) | `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` <br> `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` <br> `pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds` <br> `pbta-gene-counts-rsem-expected_count-collapsed.polya.rds` <br> `pbta-histologies.tsv`


## Color Palette Usage

This project has a set of unified color palettes.
There are 5 sets of hex color keys to be used for all final figures, stored as 5 TSV files in the `figures/palettes` folder.
`hex_codes` contains the colors to be passed to your plotting code and `color_names` contains short descriptors of each color (e.g. `gradient_1`, or `divergent_neutral`).
Each palette contains an `na_color` that is the same color in all palettes.
This color should be used for all `NA` values.
`na_color` is always the last value in the  palette.
If `na_color` is not needed or is supplied separately to a plotting function, you can use a `dplyr::filter(hex_code != "na_color")` to remove `na_color`.

To see a summary of what colors are used for histology labeling, see [`mapping-histology-labels.nb.html`](./mapping-histology-labels.nb.html)

| Palette File Name | HEX color key | Color Notes | Variable application |
|--------------|--------------------|-----------|----------------------|
|`gradient_col_palette.tsv`| <br>gradient_0:![f7f7f7](https://placehold.it/150x40/f7f7f7/000000?text=f7f7f7) <br>gradient_1:![f7fcf5](https://placehold.it/150x40/f7fcf5/000000?text=f7fcf5) <br>gradient_2:![e5f5e0](https://placehold.it/150x40/e5f5e0/000000?text=e5f5e0) <br>gradient_3:![c7e9c0](https://placehold.it/150x40/c7e9c0/000000?text=c7e9c0) <br>gradient_4:![a1d99b](https://placehold.it/150x40/a1d99b/FFFFFF?text=a1d99b) <br>gradient_5:![74c476](https://placehold.it/150x40/74c476/FFFFFF?text=74c476) <br>gradient_6:![41ab5d](https://placehold.it/150x40/41ab5d/FFFFFF?text=41ab5d) <br>gradient_7:![238b45](https://placehold.it/150x40/238b45/FFFFFF?text=238b45) <br>gradient_8:![006d2c](https://placehold.it/150x40/006d2c/FFFFFF?text=006d2c) <br>gradient_9:![00441b](https://placehold.it/150x40/00441b/FFFFFF?text=00441b) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/000000?text=f1f1f1)|10 hex_codes where gradient_0 is for an absolute `0` but may need to be removed from the palette depending on the application|For numeric data being plotted e.g. tumor mutation burden|
|`divergent_col_palette.tsv`|<br>divergent_low_5:![053061](https://placehold.it/150x40/053061/FFFFFF?text=053061) <br>divergent_low_4:![2166ac](https://placehold.it/150x40/2166ac/FFFFFF?text=2166ac) <br>divergent_low_3:![4393c3](https://placehold.it/150x40/4393c3/FFFFFF?text=4393c3) <br>divergent_low_2:![92c5de](https://placehold.it/150x40/92c5de/FFFFFF?text=92c5de) <br>divergent_low_1:![d1e5f0](https://placehold.it/150x40/d1e5f0/FFFFFF?text=d1e5f0) <br>divergent_neutral:![f7f7f7](https://placehold.it/150x40/f7f7f7/FFFFFF?text=f7f7f7) <br>divergent_high_1:![fddbc7](https://placehold.it/150x40/fddbc7/FFFFFF?text=fddbc7) <br>divergent_high_2:![f4a582](https://placehold.it/150x40/f4a582/FFFFFF?text=f4a582) <br>divergent_high_3:![d6604d](https://placehold.it/150x40/d6604d/FFFFFF?text=d6604d) <br>divergent_high_4:![b2182b](https://placehold.it/150x40/b2182b/FFFFFF?text=b2182b) <br>divergent_high_5:![67001f](https://placehold.it/150x40/67001f/FFFFFF?text=67001f) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/FFFFFF?text=f1f1f1)|12 hex codes where the numbers in the name indicate distance from `divergent_neutral`.|For data has that is bidirectional e.g. Amplification/Deletion values like `seg.mean`|
|`binary_col_palette.tsv` |<br>binary_1:![2166ac](https://placehold.it/150x40/2166ac/FFFFFF?text=2166ac) <br>binary_2:![b2182b](https://placehold.it/150x40/b2182b/FFFFFF?text=b2182b) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/000000?text=f1f1f1)|A vector of two hex codes|For binary variables e.g. presence/absence or Amp/Del as statuses|
| `oncoprint_color_palette.tsv` | <br>Missense_Mutation:![35978f](https://placehold.it/150x40/35978f/FFFFFF?text=35978f) <br>Nonsense_Mutation:![000000](https://placehold.it/150x40/000000/FFFFFF?text=000000) <br>Frame_Shift_Del:![56B4E9](https://placehold.it/150x40/56B4E9/FFFFFF?text=56B4E9) <br>Frame_Shift_Ins:![FFBBFF](https://placehold.it/150x40/FFBBFF/FFFFFF?text=FFBBFF) <br>Splice_Site:![F0E442](https://placehold.it/150x40/F0E442/FFFFFF?text=F0E442) <br>Translation_Start_Site:![191970](https://placehold.it/150x40/191970/FFFFFF?text=191970) <br>Nonstop_Mutation:![545454](https://placehold.it/150x40/545454/FFFFFF?text=545454) <br>In_Frame_Del:![CAE1FF](https://placehold.it/150x40/CAE1FF/FFFFFF?text=CAE1FF) <br>In_Frame_Ins:![FFE4E1](https://placehold.it/150x40/FFE4E1/FFFFFF?text=FFE4E1) <br>Stop_Codon_Ins:![CC79A7](https://placehold.it/150x40/CC79A7/FFFFFF?text=CC79A7) <br>Start_Codon_Del:![56B4E9](https://placehold.it/150x40/56B4E9/FFFFFF?text=56B4E9) <br>Fusion:![7B68EE](https://placehold.it/150x40/7B68EE/FFFFFF?text=7B68EE) <br>Multi_Hit:![00F021](https://placehold.it/150x40/00F021/FFFFFF?text=00F021) <br>Hom_Deletion:![313695](https://placehold.it/150x40/313695/FFFFFF?text=313695) <br>Hem_Deletion:![abd9e9](https://placehold.it/150x40/abd9e9/FFFFFF?text=abd9e9) <br>amplification:![c51b7d](https://placehold.it/150x40/c51b7d/FFFFFF?text=c51b7d) <br>loss:![0072B2](https://placehold.it/150x40/0072B2/FFFFFF?text=0072B2) <br>gain:![D55E00](https://placehold.it/150x40/D55E00/FFFFFF?text=D55E00) <br>High_Level_Gain:![FF0000](https://placehold.it/150x40/FF0000/FFFFFF?text=FF0000) <br>Multi_Hit_Fusion:![CD96CD](https://placehold.it/150x40/CD96CD/FFFFFF?text=CD96CD) | A named vector of hex codes assigned to each `short_histology` and to each `CNV`, `SNV` and `Fusion` category | For plotting an oncoprint figure, this vector provides hex codes for `CNV`, `SNV`, and `Fusion` categories |
| `tumor_descriptor_palette.tsv` | <br> Initial CNS Tumor:	![709AE1FF](https://placehold.it/150x40/35978f/FFFFFF?text=709AE1FF)  <br> Progressive	![075149FF](https://placehold.it/150x40/35978f/FFFFFF?text=709AE1FF) <br> Progressive Disease Post-Mortem	![075149FF](https://placehold.it/150x40/35978f/FFFFFF?text=075149FF)
<br> Recurrence	![FD8CC1FF](https://placehold.it/150x40/35978f/FFFFFF?text=FD8CC1FF) <br> Second Malignancy	![FD7446FF](https://placehold.it/150x40/35978f/FFFFFF?text=FD7446FF) | A named vector of hex codes assigned to each `tumor_descriptor` | For plotting in sample distribution, this vector provides color for tumor descriptor categories |

## Color coding examples in R

#### Example 1) Color coding by histologies

**Step 1)** Read in color palette and select the pertinent columns

There's some extra columns in `histology_label_color_table.tsv` that you don't need for plotting per se but are more record-keeping purposes. 
With the code chunk below, we only import the four columns we need and then do a factor reorder to make sure the `display_group` is in the order declared by `display_order`. 

```
# Import standard color palettes for project
histology_label_mapping <- readr::read_tsv(
  file.path(figures_dir, "palettes", "histology_label_color_table.tsv")
  ) %>% 
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>% 
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))
```

**Step 2)** Use `dplyr::inner_join` using `Kids_First_Biospecimen_ID` to join by so you can add on the `hex_codes` and `display_group` for each biospecimen. 
`display_order` specifies what order the `display_group`s should be displayed.

```
# Read in the metadata
metadata <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID")
```

**Step 4)** Make your plot and use the `hex_codes` and `display_group` columns.  

Using the `ggplot2::scale_fill_identity()` or `ggplot2::scale_color_identity()` allows you to supply the `hex_codes` column from a color palette to `ggplot2` with a `fill` or `color` argument respectively.
For base R plots, you should be able to supply the `hex_codes` column as your `col` argument.
`display_group` should be used as the labels in the plot.

```
metadata %>%
  dplyr::group_by(display_group, hex_codes) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  ggplot2::ggplot(ggplot2::aes(x = display_group, y = count, fill = hex_codes)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_identity()
```

#### Example 2) Color coding by numeric data

**Step 1)** Import the palette.

You may want to remove the `na_color` at the end of the list depending on whether your data include `NA`s or if the plotting function you are using has the `na_color` supplied separately.

```
gradient_col_palette <- readr::read_tsv(
  file.path(figures_dir, "palettes", "gradient_color_palette.tsv")
)
```

If we need the `NA` color separated, like for use with `ComplexHeatmap` which has a separate argument for the color for `NA` values.

```
na_color <- gradient_col_palette %>%
  dplyr::filter(color_names == "na_color")

gradient_col_palette <- gradient_col_palette %>%
  dplyr::filter(color_names != "na_color")
```

**Step 2)** Make a color function.  

In this example, we are building a `colorRamp2` function based on a regular interval between the minimum and maximum of our variable `df$variable` by using `seq`.
However, depending on your data's distribution a regular interval based palette might not represent your data well on the plot.
You can provide any numeric vector to color code a palette using `circlize::colorRamp2` as long as that numeric vector is the same length as the palette itself.

```
gradient_col_val <- seq(from = min(df$variable), to = max(df$variable),
                        length.out = nrow(gradient_col_palette))

col_fun <- circlize::colorRamp2(gradient_col_val,
                                gradient_col_palette$hex_codes)
```
**Step 3)** Apply to numeric data, or supply to your plotting code.  

This step depends on how your main plotting function would like the data supplied.
For example, `ComplexHeatmap` wants a function to be supplied to their `col` argument.

```
# Apply to variable directly and make a new column
df <- df %>%
  dplyr::mutate(color_key = col_fun(variable))

## OR ##

# Some plotting packages want a color function

ComplexHeatmap::Heatmap(
  df,
  col = col_fun,
  na_col = na_color$hex_codes
)
```

### Updating color palettes

The non-histologies color palette TSV files are created by running `scripts/color_palettes.R`, which can be called by `Rscript scripts/color_palettes.R`.
Hex codes for the palettes are hard-coded in this script.
The script can be called from anywhere in this repository (will look for the `.git` file).
The hex codes table in `figures/README.md` and its swatches should also be updated by using the `swatches_table` function at the end of the script and copy and pasting this function's output to the appropriate place in the table.

The histology color palette file is created by running `Rscript -e "rmarkdown::render('figures/mapping-histology-labels.Rmd', clean = TRUE)"`.
