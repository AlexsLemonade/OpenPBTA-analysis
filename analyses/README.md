## Analysis Modules

This directory contains various analysis modules in the OpenPBTA project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance

The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
In addition, this table reflects which analyses are included in the OpenPBTA manuscript.
This is in service of documenting interdependent analyses.
Note that _nearly all_ modules use the harmonized clinical data file (`pbta-histologies.tsv`) even when it is not explicitly included in the table below.

| [`survival-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis) | `pbta-histologies.tsv` <br> `independent-specimens.wgswxs.primary.tsv` (from results by independent-samples module) <br> `TelomeraseScores_PTBAStranded_FPKM.txt` (from results by telomerase-activity-prediction module) <br> `tp53_altered_status.tsv` (from results by tp53_nf1_score module) |  `cox_reg_results_per_{covariates}.tsv` (see module results for details) | N/A
