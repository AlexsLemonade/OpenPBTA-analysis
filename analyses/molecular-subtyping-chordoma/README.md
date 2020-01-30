## Molecular subtyping of chordomas

**Module authors:** Mateusz Koptyra ([@mkoptyra](https://github.com/mkoptyra))

This module consists of a single notebook that looks at _SMARCB1_ focal copy status and expression levels.
It can be run via the command line with the following:

```
Rscript -e "rmarkdown::render('01-Subtype-chordoma.Rmd', clean = TRUE)"
```

### Notes on copy status

This notebook uses an older version of annotated CNVkit file from the `focal-cn-file-preparation` ([`fa21429`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/fa214291713575be7fd20c92374b268870f4173f)) as the current version of the annotated file from CNVkit may be too restrictive (see: [#473](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/473)). 
This will need to be updated to use the copy number consensus file as well.
