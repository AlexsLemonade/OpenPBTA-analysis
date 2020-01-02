## Plotting GISTIC results

**Module Author:** Candace Savonen ([@cansavvy](https://www.github.com/cansavvy))

The goal of this analysis is to plot GISTIC results and make CNV plots by histology groups. 

### Running the analysis

This analysis consists of a single R Notebook, that can be run with the following from the top directory of the project:

```
Rscript -e "rmarkdown::render('analyses/cnv-chrom-plot/gistic_plot.Rmd', clean = TRUE)"
```

### Output

The output is a plot of the GISTIC scores (`plots/gistic_plot.png`) as well as
plots of the `seg.mean` by each histology group (e.g. `plots/Chondrosarcoma_plot.png`).
