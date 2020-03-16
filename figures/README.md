## Color Palette Usage

This project has a set of unified color palettes.
There are six sets of hex color keys to be used for all final figures.
They are created by running `scripts/color-palettes.R` which creates 6 TSV files
in the `figures/palettes` folder.
Each palette TSV file has two columns of data: `hex_codes` and `color_names`.
`hex_codes` contains the colors to be passed to your plotting code and `color_names` contains short descriptors of each color (e.g. `gradient_1`, or `divergent_neutral`).
Each palette contains an `na_color` that is the same color in all palettes it
is always the last value in the list.
So depending on the circumstance, if `na_color` is not needed or will be supplied
to a plot's argument in a different way, you can use a `dplyr::filter(hex_code != "na_color")`
to remove the `na_color`.
Biospecimens without a `short_histology` designation are assigned the `NA` color
for this project.

| Palette File Name | HEX color key | Color Notes | Variable application |
|--------------|--------------------|-----------|----------------------|
|`histology_color_palette.tsv`|<br>Adenoma:![f23d3d](https://placehold.it/150x40/f23d3d/FFFFFF?text=f23d3d) <br>ATRT:![731d1d](https://placehold.it/150x40/731d1d/FFFFFF?text=731d1d) <br>Central neurocytoma:![b38686](https://placehold.it/150x40/b38686/FFFFFF?text=b38686) <br>Chondrosarcoma:![cc5c33](https://placehold.it/150x40/cc5c33/FFFFFF?text=cc5c33) <br>Chordoma:![331c0d](https://placehold.it/150x40/331c0d/FFFFFF?text=331c0d) <br>Choroid plexus tumor:![ffb380](https://placehold.it/150x40/ffb380/FFFFFF?text=ffb380) <br>CNS EFT-CIC:![b25f00](https://placehold.it/150x40/b25f00/FFFFFF?text=b25f00) <br>CNS lymphoma:![f2d6b6](https://placehold.it/150x40/f2d6b6/000000?text=f2d6b6) <br>CNS neuroblastoma:![736556](https://placehold.it/150x40/736556/FFFFFF?text=736556) <br>CNS Rhabdomyosarcoma:![ffaa00](https://placehold.it/150x40/ffaa00/FFFFFF?text=ffaa00) <br>CNS sarcoma:![4c3d00](https://placehold.it/150x40/4c3d00/FFFFFF?text=4c3d00) <br>Craniopharyngioma:![e2f200](https://placehold.it/150x40/e2f200/FFFFFF?text=e2f200) <br>DNET:![919926](https://placehold.it/150x40/919926/FFFFFF?text=919926) <br>Dysplasia:![d6f2b6](https://placehold.it/150x40/d6f2b6/000000?text=d6f2b6) <br>Embryonal Tumor:![304d26](https://placehold.it/150x40/304d26/FFFFFF?text=304d26) <br>Ependymoma:![00f241](https://placehold.it/150x40/00f241/FFFFFF?text=00f241) <br>ETMR:![009929](https://placehold.it/150x40/009929/FFFFFF?text=009929) <br>Ganglioglioma:![698c7c](https://placehold.it/150x40/698c7c/FFFFFF?text=698c7c) <br>Germinoma:![39e6c3](https://placehold.it/150x40/39e6c3/FFFFFF?text=39e6c3) <br>Glial-neuronal tumor NOS:![005359](https://placehold.it/150x40/005359/FFFFFF?text=005359) <br>Gliosis:![263233](https://placehold.it/150x40/263233/FFFFFF?text=263233) <br>Hemangioblastoma:![00c2f2](https://placehold.it/150x40/00c2f2/FFFFFF?text=00c2f2) <br>Hemangioma:![40a6ff](https://placehold.it/150x40/40a6ff/FFFFFF?text=40a6ff) <br>HGAT:![406280](https://placehold.it/150x40/406280/FFFFFF?text=406280) <br>Langerhans Cell histiocytosis:![0044ff](https://placehold.it/150x40/0044ff/FFFFFF?text=0044ff) <br>LGAT:![00144d](https://placehold.it/150x40/00144d/FFFFFF?text=00144d) <br>LGMT:![acbbe6](https://placehold.it/150x40/acbbe6/FFFFFF?text=acbbe6) <br>Medulloblastoma:![7373e6](https://placehold.it/150x40/7373e6/FFFFFF?text=7373e6) <br>Meningioma:![3d0099](https://placehold.it/150x40/3d0099/FFFFFF?text=3d0099) <br>MPNST:![c200f2](https://placehold.it/150x40/c200f2/FFFFFF?text=c200f2) <br>Neurofibroma:![917399](https://placehold.it/150x40/917399/FFFFFF?text=917399) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/000000?text=f1f1f1) <br>Oligodendroglioma:![f279da](https://placehold.it/150x40/f279da/FFFFFF?text=f279da) <br>Other:![cc0052](https://placehold.it/150x40/cc0052/FFFFFF?text=cc0052) <br>Pineoblastoma:![994d6b](https://placehold.it/150x40/994d6b/FFFFFF?text=994d6b) <br>Schwannoma:![4d2636](https://placehold.it/150x40/4d2636/FFFFFF?text=4d2636) <br>Teratoma:![ffbfd9](https://placehold.it/150x40/ffbfd9/FFFFFF?text=ffbfd9)|a named vector of the hex values that were assigned to each `short_histology` group table|For color-coding by `short_histology` when it's more convenient to assign colors by `short_histology` category.|
|`gradient_col_palette.tsv`| <br>gradient_0:![f7f7f7](https://placehold.it/150x40/f7f7f7/000000?text=f7f7f7) <br>gradient_1:![f7fcf5](https://placehold.it/150x40/f7fcf5/000000?text=f7fcf5) <br>gradient_2:![e5f5e0](https://placehold.it/150x40/e5f5e0/000000?text=e5f5e0) <br>gradient_3:![c7e9c0](https://placehold.it/150x40/c7e9c0/000000?text=c7e9c0) <br>gradient_4:![a1d99b](https://placehold.it/150x40/a1d99b/FFFFFF?text=a1d99b) <br>gradient_5:![74c476](https://placehold.it/150x40/74c476/FFFFFF?text=74c476) <br>gradient_6:![41ab5d](https://placehold.it/150x40/41ab5d/FFFFFF?text=41ab5d) <br>gradient_7:![238b45](https://placehold.it/150x40/238b45/FFFFFF?text=238b45) <br>gradient_8:![006d2c](https://placehold.it/150x40/006d2c/FFFFFF?text=006d2c) <br>gradient_9:![00441b](https://placehold.it/150x40/00441b/FFFFFF?text=00441b) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/000000?text=f1f1f1)|10 hex_codes where gradient_0 is for an absolute `0` but may need to be removed from the palette depending on the application|For numeric data being plotted e.g. tumor mutation burden|
|`divergent_col_palette.tsv`|<br>divergent_low_5:![053061](https://placehold.it/150x40/053061/FFFFFF?text=053061) <br>divergent_low_4:![2166ac](https://placehold.it/150x40/2166ac/FFFFFF?text=2166ac) <br>divergent_low_3:![4393c3](https://placehold.it/150x40/4393c3/FFFFFF?text=4393c3) <br>divergent_low_2:![92c5de](https://placehold.it/150x40/92c5de/FFFFFF?text=92c5de) <br>divergent_low_1:![d1e5f0](https://placehold.it/150x40/d1e5f0/FFFFFF?text=d1e5f0) <br>divergent_neutral:![f7f7f7](https://placehold.it/150x40/f7f7f7/FFFFFF?text=f7f7f7) <br>divergent_high_1:![fddbc7](https://placehold.it/150x40/fddbc7/FFFFFF?text=fddbc7) <br>divergent_high_2:![f4a582](https://placehold.it/150x40/f4a582/FFFFFF?text=f4a582) <br>divergent_high_3:![d6604d](https://placehold.it/150x40/d6604d/FFFFFF?text=d6604d) <br>divergent_high_4:![b2182b](https://placehold.it/150x40/b2182b/FFFFFF?text=b2182b) <br>divergent_high_5:![67001f](https://placehold.it/150x40/67001f/FFFFFF?text=67001f) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/FFFFFF?text=f1f1f1)|12 hex codes where the numbers in the name indicate distance from `divergent_neutral`.|For data has that is bidirectional e.g. Amplification/Deletion values like `seg.mean`|
|`binary_col_palette.tsv` |<br>binary_1:![2166ac](https://placehold.it/150x40/2166ac/FFFFFF?text=2166ac) <br>binary_2:![b2182b](https://placehold.it/150x40/b2182b/FFFFFF?text=b2182b) <br>na_color:![f1f1f1](https://placehold.it/150x40/f1f1f1/000000?text=f1f1f1)|A vector of two hex codes|For binary variables e.g. presence/absence or Amp/Del as statuses|


## Color coding examples in R

#### Example 1) Color coding by `short_histology`.

**Step 1)** Read in color palette and format as a named list
```
histology_col_palette <- readr::read_tsv(
  file.path("figures", "palettes", "histology_color_palette.tsv")
  ) %>%
  # We'll use deframe so we can use it as a recoding list
  tibble::deframe()
```

**Step 2)** For any data.frame with a `short_histology` column, recode NAs as "none".
```
metadata <- readr::read_tsv(file.path("data", "pbta-histologies.tsv") %>%
  # Easier to deal with NA short histologies if they are labeled something different
  dplyr::mutate(short_histology = as.character(tidyr::replace_na(short_histology, "none")))
```

**Step 3)** Use dplyr::recode on `short_histology` column to make a new color column.
```
metadata <- metadata %>%
  # Tack on the sample color using the short_histology column and a recode
  dplyr::mutate(sample_color = dplyr::recode(short_histology,
                                             !!!histology_col_palette))
```

**Step 4)** Make your plot and use the `sample_color` column.  
Using the `ggplot2::scale_fill_identity()` or `ggplot2::scale_color_identity()`
allows you to supply an exact `hex_code` column to `ggplot2` with a `fill` or
`color` argument respectively.
For base R plots, you should be able to supply the `sample_color` column as your
`col` argument.
```
metadata %>%
  dplyr::group_by(short_histology, sample_color) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  ggplot2::ggplot(ggplot2::aes(x = short_histology, y = count, fill = sample_color)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_identity()
```

#### Example 2) Color coding by numeric data

**Step 1)** Import the palette.
```
gradient_col_palette <- readr::read_tsv(
  file.path(figures_dir, "palettes", "gradient_color_palette.tsv")
  ) %>%
  # We won't need NA color in this instance, ComplexHeatmap has a separate argument for that
  dplyr::filter(color_names != "na_color")
```
**Step 2)** Make a color function.

Note that the last value in every data.frame is the `NA` color.
The numbers supplied here will be highly dependent on what your data's distribution looks like.
```
gradient_col_val <- seq(from = min(df$variable), to = max(df$variable),
                        length.out = length(gradient_col_palette))

col_fun <- circlize::colorRamp2(gradient_col_val,
                                gradient_col_palette)
```
**Step 3)** Apply to numeric data, or supply to your plotting code
This step depends on how your main plotting function would like the data supplied.
For example, ComplexHeatmap wants a user to supply a `function`.
```
# Apply to variable directly and make a new column
df <- df %>%
  dplyr::mutate(color_key = col_fun(variable))

## OR ##

# Some plotting packages want a color function
ComplexHeatmap::heatmap(df,
  col = col_fun)
```
