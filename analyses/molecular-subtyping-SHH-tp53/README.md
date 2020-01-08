## _TP53_ mutation status in Medulloblastoma SHH subtype samples

**Module Author:** Candace Savonen ([@cansavvy](https://www.github.com/cansavvy))

The goal of this analysis is to further classify SHH subtype medulloblastoma samples into SHH, _TP53_-mutated and SHH, _TP53_-wildtype per [#247](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/247).

We use `molecular_subtype` information from the harmonized clinical file (`pbta-histologies.tsv`) to restrict our analysis to samples classified as `SHH`.
The medulloblastoma subtype classifier uses RNA-seq data, so we must identify WGS biospecimens that map to the same `sample_id`, if available.
We then look at the presence or absence of mutations in _TP53_ in the consensus mutation file (`pbta-snv-consensus-mutation.maf.tsv.gz`).

### Running the analysis

This analysis consists of a single R Notebook, that can be run with the following from the top directory of the project:

```
Rscript -e "rmarkdown::render('analyses/molecular-subtyping-SHH-tp53/SHH-tp53-molecular-subtyping-data-prep.Rmd', clean = TRUE)"
```

### Output

The output is a table _TP53_ mutation status on a per `sample_id` basis available as `results/tp53-shh-samples-status.tsv`.
