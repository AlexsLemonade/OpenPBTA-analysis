# Figures

This directory contains the code required to generate individual panels for main and supplemental display items.
It is designed such that a single shell script will run all required data processing or analysis steps prior to plotting.
There are two main strategies for creating individual panels: 1. scripts that are exclusively created for publication-ready display (located in `figures/scripts`) or 2. copying publication-ready plots that are generated within a given analysis module (i.e., `analyses`) to the correct location in `figures/`.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

  - [Running the figure generation script](#running-the-figure-generation-script)
      - [1. Obtain the current dataset.](#1-obtain-the-current-dataset)
      - [2. Set up an up-to-date project Docker container.](#2-set-up-an-up-to-date-project-docker-container)
      - [3. Run the bash script that generates the figures (`scripts/generate-figures.sh`).](#3-run-the-bash-script-that-generates-the-figures-scriptsgenerate-figuressh)
  - [Adding or updating figures](#adding-or-updating-figures)
    - [Individual panels and where to find them](#individual-panels-and-where-to-find-them)
    - [Compiled multipanel figures](#compiled-multipanel-figures)
    - [Documenting individual figures & scripts](#documenting-individual-figures--scripts)
  - [Figure Guidelines](#figure-guidelines)
    - [Color Palettes](#color-palettes)
      - [Updating color palettes](#updating-color-palettes)
      - [Examples palette usage in R](#examples-palette-usage-in-r)
        - [Example 1) Color coding by disease label in `ggplot2`](#example-1-color-coding-by-disease-label-in-ggplot2)
        - [Example 2) Color coding by numeric data](#example-2-color-coding-by-numeric-data)
    - [Overall figure theme](#overall-figure-theme)
    - [Statistics](#statistics)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running the figure generation script

Follow these steps to refresh all publication-ready figures.
All steps assume your current directory is the top of this repository.

#### 1. Obtain the current dataset.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#how-to-obtain-openpbta-data) to obtain the current data release.
We recommend [using the download script](https://github.com/AlexsLemonade/OpenPBTA-analysis#data-access-via-download-script) to obtain data because this will automatically create symlinks in `data/` to the latest files.

#### 2. Set up an up-to-date project Docker container.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-image) for setting up the project Docker container.
Briefly, the latest version of the project Docker image, which is updated upon commit to `master`, can be obtained and run via:

```bash
docker pull ccdlopenpbta/open-pbta:latest
docker run \
  -e PASSWORD=<password> \
  -p 8787:8787 \
  -v $(pwd):/home/rstudio/kitematic \
  ccdlopenpbta/open-pbta:latest
```
You may choose to use [`docker exec`](https://docs.docker.com/engine/reference/commandline/exec/) to interact with the container from there or if you'd prefer the RStudio interface, you can navigate to `localhost:8787` and enter username `rstudio` and the password you set with the `run` command above.

#### 3. Run the bash script that generates the figures (`scripts/generate-figures.sh`).

_⚠️ This requires 64GB of RAM to run successfully! You can set `RUN_LOCAL=1` for local testing that skips the RAM-intensive steps._

This script runs **_all_** the intermediate steps needed to generate figures starting with the files in the data release.
Earlier versions of this script ran the consensus SNV modules, but the consensus SNV files are included in the data download and the downstream modules use the files in the `data/` directory, so this is no longer necessary.

```bash
bash figures/generate-figures.sh
```

In some cases, the shell script will copy plots to the appropriate locations in the `figures/` directory.

## Adding or updating figures

All figures need to adhere to [the guidelines in this README](#figure-guidelines) and need to be generated, including any upstream steps, via `figures/generate-figures.sh`.
This section includes details for file organization and guidelines that are not covered in the section on style.

### Individual panels and where to find them

Individual panels, including any standalone legends, should be saved as **PDFs** in `figures/pdfs/<figure>/panels`.
Each figure has its own subdirectory within `figures/pdfs`; individual panels for Figure 1 would be saved in `figures/pdfs/fig1/panels`, for Figure 2 in `figures/pdfs/fig2/panels` and so on and so forth.  

### Compiled multipanel figures

Multipanel figures are compiled with Adobe Illustrator and saved as PDFs in their respective subdirectories.
For example, the compiled, PDF version of Figure 4 is saved as `figures/pdf/fig4/figure4.pdf`.
(Original `.ai` files are available in Google Drive and can be shared as necessary -- contact @jaclyn-taroni!)

PNG versions of figures are exported using the PDF file and saved to the `figures/pngs` folder.
PNGs version can be linked in the accompanying manuscript repository [`AlexsLemonade/OpenPBTA-manuscript`](https://github.com/AlexsLemonade/OpenPBTA-manuscript/).

### Documenting individual figures & scripts

Each figure should include its own associated README in the appropriate PDF subdirectory, e.g., the documentation for Figure 4 can be found in `figures/pdfs/fig4/README.md`.
That README should detail what analysis modules or individual scripts are required to generate the individual panels for that figure.

Scripts that are exclusively for creating publication-ready plots are available in `figures/scripts` and should be documented in the `figures/scripts/README.md`.

## Figure Guidelines

### Color Palettes 

This project has a set of unified color palettes.
There are 6 sets of hex color keys to be used for all final figures, stored as 6 TSV files in the `figures/palettes` folder.
`hex_codes` contains the colors to be passed to your plotting code and `color_names` contains short descriptors of each color (e.g. `gradient_1`, or `divergent_neutral`).

Palettes for numeric or binary data contain an `na_color` that is uniform across palettes.
This color should be used for all `NA` values.
`na_color` is always the last value in the  palette.
If `na_color` is not needed or is supplied separately to a plotting function, you can use a `dplyr::filter(hex_code != "na_color")` to remove `na_color`.

To see a summary of what colors are used for histology labeling, see [`mapping-histology-labels.nb.html`](./mapping-histology-labels.nb.html)

| Palette File Name | Color Notes | Variable application |
|--------------|-----------|----------------------|
|`gradient_col_palette.tsv` | 10 hex_codes where gradient_0 is for an absolute `0` but may need to be removed from the palette depending on the application| For numeric data being plotted e.g., tumor mutation burden |
|`divergent_col_palette.tsv` | 12 hex codes where the numbers in the name indicate distance from `divergent_neutral`. | For data has that is bidirectional e.g., Amplification/Deletion values like `seg.mean`|
|`binary_col_palette.tsv` | A vector of two hex codes | For binary variables e.g., presence/absence or Amp/Del as statuses |
| `oncoprint_color_palette.tsv` | A named vector of hex codes assigned to each `short_histology` and to each `CNV`, `SNV` and `Fusion` category | For plotting an oncoprint figure, this vector provides hex codes for specific alterations |
| `tumor_descriptor_palette.tsv` | A named vector of hex codes assigned to each `tumor_descriptor` | For plotting sample distributions, this vector provides color for tumor descriptor categories |
|`broad_histology_cancer_group_palette.tsv` | Contains multiple columns having to do with the display by disease label (i.e., `broad_histology` or `cancer_group`) | To be used for any plots that require coloring by `broad_histology` or `cancer_group`; please see `figures/mapping-histology-labels.Rmd` for more information |

_⚠️ `histology_label_color_table.tsv` is a deprecated version of the palettes used for `broad_histology` and `cancer_group`._ 
_It has been retained, for the moment, to prevent the introduction of breaking changes to multiple analysis modules._
_It is not to be used for any future development and the relationships between individual biospecimens and disease labels may be out of date._

#### Updating color palettes

The non-histologies color palette TSV files are created by running `scripts/color_palettes.R`, which can be called by `Rscript scripts/color_palettes.R`.
Hex codes for the palettes are hard-coded in this script.
Do not manually edit palette TSV files.
The script can be called from anywhere in this repository (will look for the `.git` file).
The hex codes table in `figures/README.md` and its swatches should also be updated by using the `swatches_table` function at the end of the script and copy and pasting this function's output to the appropriate place in the table.

The histology color palette file is created by running `Rscript -e "rmarkdown::render('figures/mapping-histology-labels.Rmd', clean = TRUE)"`.

#### Examples palette usage in R

##### Example 1) Color coding by disease label in `ggplot2`

**Step 1)** Read in color palette file and create a named vector of hex codes

Here's an example for `cancer_group`, specifically:

```r
# Get palette for cancer group
cancer_group_palette <- readr::read_tsv(
  "figures/palettes/broad_histology_cancer_group_palette.tsv")
) %>%
  dplyr::select(cancer_group, cancer_group_hex) %>%
  # Remove NA values -- a cancer group hex value will be NA only if the
  # cancer group is NA
  dplyr::filter(complete.cases(.))

# Make color palette suitable for use with ggplot
annotation_colors <- cancer_group_palette$cancer_group_hex
names(annotation_colors) <- cancer_group_palette$cancer_group  
```

**Step 2)** Make your ggplot using the named vector as a manual palette

You will be able to use the named vector with `ggplot2` functions such as `scale_fill_manual()` or `scale_color_manual()`, like so:

```r
ggplot2::ggplot(
  ggplot2::aes(x = cancer_group,
               y = y_value, 
               fill = cancer_group)
) +
  geom_boxplot() +
  scale_fill_manual(values = annotation_colors)
```

##### Example 2) Color coding by numeric data

**Step 1)** Import the palette.

You may want to remove the `na_color` at the end of the list depending on whether your data include `NA`s or if the plotting function you are using has the `na_color` supplied separately.

```r
gradient_col_palette <- readr::read_tsv(
  file.path(figures_dir, "palettes", "gradient_color_palette.tsv")
)
```

If we need the `NA` color separated, like for use with `ComplexHeatmap` which has a separate argument for the color for `NA` values.

```r
na_color <- gradient_col_palette %>%
  dplyr::filter(color_names == "na_color")

gradient_col_palette <- gradient_col_palette %>%
  dplyr::filter(color_names != "na_color")
```

**Step 2)** Make a color function.  

In this example, we are building a `colorRamp2` function based on a regular interval between the minimum and maximum of our variable `df$variable` by using `seq`.
However, depending on your data's distribution a regular interval based palette might not represent your data well on the plot.
You can provide any numeric vector to color code a palette using `circlize::colorRamp2` as long as that numeric vector is the same length as the palette itself.

```r
gradient_col_val <- seq(from = min(df$variable), to = max(df$variable),
                        length.out = nrow(gradient_col_palette))

col_fun <- circlize::colorRamp2(gradient_col_val,
                                gradient_col_palette$hex_codes)
```

**Step 3)** Apply to numeric data, or supply to your plotting code.  

This step depends on how your main plotting function would like the data supplied.
For example, `ComplexHeatmap` wants a function to be supplied to their `col` argument.

```r
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

### Overall figure theme

In general, we will use the `ggpubr` package with `ggtheme = theme_pubr())` and color palette `simpsons` from package `ggsci` since it has 16 levels and can accommodate the levels in groups such as `molecular_subtype`.

To view the palette:

```r
scales::show_col(ggsci::pal_simpsons("springfield")(16))
```

For 2+ group comparisons, we will use violin or boxplots with jitter.

### Statistics

Some modules perform group-wise comparisons. 
For the manuscript, we may want to output tables of the statistics and/or print the statistical test and p-value directly on the plot.
We use the functions `ggpubr::compare_means()` and `ggpubr::stat_compare_means()` for this. 
Below are the default tests, parameters, and method options for 2 groups or [more than two groups](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/#compare-more-than-two-groups) for your convenience.
Caution: the default p-values on the plots are uncorrected.

|                                            | 2 groups                                             | 3+ groups                                                             |
|--------------------------------------------|------------------------------------------------------|-----------------------------------------------------------------------|
| Default test (method)                      | Wilcoxon                                             | Kruskal-wallis                                                        |
| Allowed methods                            | "wilcox.test" (non-parametric) "t.test" (parametric) | "kruskal.test" (non-parametric) "anova" (parametric)                  |
| Default multiple testing (p.adjust.method) | NA                                                   | yes, but not bonferroni                                               |
| Allowed p.adjust.method                    | NA                                                   | "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" |

Below is an example for creating a violin plot with boxplot, jitter, and appropriate statistics.

```r
if(length(unique(df$var_x)) > 2){
    method <- "kruskal.test"
  } else {
    method <- "wilcox.test"
  }


p <- ggviolin(df, x = "var_x", y = "var_y", 
           color = "var_color", 
           palette = "simpsons",
           order = c("a", "b", "c"),
           add = c("boxplot", "jitter"),  
           ggtheme = theme_pubr()) +
    # Add pairwise comparisons p-value
    stat_compare_means(method = method, label.y = 1.2, label.x.npc = "center") +
    xlab("xlab_text") +
    ylab("ylab_text") +
    rremove("legend")
```
