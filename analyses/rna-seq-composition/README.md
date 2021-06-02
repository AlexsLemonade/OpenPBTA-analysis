# RNA-Seq sample composition analysis

**Module authors:** Holly Beale ([@hollybeale](https://github.com/hbeale))

## Usage


## Folder content

This folder contains the script that analyzes MEND output.

The results directory contains the file samples_with_flags.tsv which contains all samples with a sample compsotions outside of a reference range or with fewer than 10 million MEND reads. 


## Usage

To re-run this analysis from the command line, run the following:
```

Rscript -e "rmarkdown::render('analyses/rna-seq-composition/rna-seq-composition.Rmd', clean = TRUE)"

```


