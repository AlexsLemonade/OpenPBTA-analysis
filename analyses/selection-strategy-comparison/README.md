## Explore RNA-seq selection strategies

**Module author:** Joshua Shapiro ([@jashapiro](https://github.com/jashapiro))

This analysis was originally undertaken prior to [`release-v5-20190924`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/release-notes.md#release-v5-20190924). 
Data releases prior to v5 had samples from the two library strategies (stranded, poly-A) in the same file, but due to the findings of this analysis, libraries of different strategies were split into separate files in the subsequent data releases. 
The notebook plots UMAP results following various normalization strategies, which can be found in `plots`.

#### Running the analysis

This notebook can be run with the following (from the root directory of the repository):

```
Rscript -e "rmarkdown::render('analyses/selection-strategy-comparison/01-selection-strategies.rmd', clean = TRUE)"
```
