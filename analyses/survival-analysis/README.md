## Perform survival analysis for samples with 

**Module authors:** Run Jin [@runjin326](https://github.com/runjin326), Stephanie J. Spielman, Jo Lynne Rokita 

In this module, we perform survival analyses by Kaplan-Meier, log-rank, and/or or cox regression for samples in OpenPBTA using different covariates. 
Detailed information about covariates can be found below.


### Usage
```sh
bash run-survival.sh
```

### Module contents

#### Function
`util/survival_models.R` performs kaplan-meier, log rank, or cox regression for univariate and/or multivariate analyses.
This function returns a list with three objects:
1. The original modeel fit object
2. The summary table
3. The original dataframe

#### Scripts

`survival-analysis_template.Rmd` showed the usage of function `survival_analysis` to perform Kaplan-Meier, log-rank and cox regression survival analyses. 
Particularly, it demonstrated runninng Kaplan-Meier and log-rank survival analysis on categorial variable, using `germline_sex_estimate` as an example.
It also demonstrated how to run cox regression survival analysis on continuous variable, using `tmb` (tumor mutation burden) as an example.

`survival-analysis_subtypes.Rmd` runs the following analyses:
1. Univariate models
- Kaplan-Meier survival analysis on HGG samples (minus oligodendrogliomas) by `molecular_subtype`
- Cox regression univariate analysis on HGG samples (minus oligodendrogliomas) by `molecular_subtype`

`survival-analysis_tp53_telomerase.Rmd` runs the following analysis:
1. Univariate models 
- extent of tumor resection (Total, Partial, Biopsy, Unknown)
- LGG group (LGG vs. non-LGG)
- HGG group (HGG vs. non-HGG)
- TP53 classifier score (as a continuous variable)
- telomerase score (EXTEND, as a continuous variable)
- Cancer group (cox)
- Cancer group (log rank)

2. Multivariate models
- Interaction model: `TP53 classifier score * telomerase score * extent of tumor resection * LGG group`
- Additive model: `TP53 classifier score + telomerase score + extent of tumor resection + LGG group`
- Interaction model: `TP53 classifier score * telomerase score * extent of tumor resection * HGG group`
- Additive model: `TP53 classifier score + telomerase score + extent of tumor resection + HGG group`
- Additive models: `TP53 classifier score + telomerase score + extent of tumor resection` for each cancer group

3. Fit results are saved as `.RDS` files for each analysis

`survival-analysis_immune.Rmd` runs the following analyses:
1. Univariate models:
- PD-L1 expression (stranded samples only)

2. Multivariate models:
- CD274 (PD-L1) expression, controlling for RNA library
- quanTIseq immune cell fractions + extent of tumor resection + LGG group across entire PBTA cohort
- quanTIseq immune cell fractions + CD274 expression + extent of tumor resection + LGG group across entire PBTA cohort
- quanTIseq immune cell fractions + CD274 expression + extent of tumor resection across suitable cancer groups (DECEASED N>=3, thus no LGG)


#### Output

```
results/*.rds
results/*tsv
plots/KM_hgg_subtypes.pdf
```

The RDS and TSV files contain survival model fit results.
The Kaplan-meier curve for the log rank analysis of HGG by molecular subtype is saved as `KM_hgg_subtypes.pdf`
