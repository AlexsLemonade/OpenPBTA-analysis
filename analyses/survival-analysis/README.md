## Perform survival analysis for samples with 

**Module authors:** Candace Savonen (ALSF CCDL), Run Jin (D3b), Stephanie J. Spielman (ALSF CCDL), Jo Lynne Rokita (D3b)

In this module, we perform survival analyses by Kaplan-Meier, log-rank, and/or cox regression for samples in OpenPBTA using different covariates. 
Detailed information about covariates can be found below.


### Usage
```sh
bash run-survival.sh
```

### Module contents

#### Function
`util/survival_models.R` performs kaplan-meier, log rank, or cox regression for univariate and/or multivariate analyses.
- The `survival_analysis` function performs survival modeling and returns a list with three objects:
    1. The original model fit object
    2. The summary table
    3. The original dataframe
- The `fit_save_model` function wraps the function `survival_analysis()` to perform and export a survival model


#### Scripts

`survival-analysis_template.Rmd` shows the usage of function `survival_analysis` to perform Kaplan-Meier, log-rank and cox regression survival analyses. 
Particularly, it demonstrates runninng Kaplan-Meier and log-rank survival analysis on a categorical variable, using `germline_sex_estimate` as an example.
It also demonstrates how to run cox regression survival analysis on a continuous variable, using `tmb` (tumor mutation burden) as an example.

`survival-analysis_subtypes.Rmd` runs the following analyses:
1. Univariate models
- Kaplan-Meier survival analysis on HGG samples by `molecular_subtype`
- Cox regression univariate analysis on HGG samples by `molecular_subtype`

`survival-analysis_tp53_telomerase.Rmd` runs the following analysis:
1. Univariate models 
- extent of tumor resection (cox)
- LGG broad histology (cox)
- HGG broad histology (cox)
- Cancer group (cox)
- Cancer group (log rank)

2. Multivariate models
- Interaction model: `TP53 classifier score * telomerase score * extent of tumor resection * LGG group` across entire PBTA cohort
- Additive model: `TP53 classifier score + telomerase score + extent of tumor resection + LGG group + HGG_group` across entire PBTA cohort
- Additive models: `TP53 classifier score + telomerase score + extent of tumor resection` for each `cancer_group`

3. Fit results are saved as `.RDS` files for each analysis

`survival-analysis_immune.Rmd` runs the following analyses:
1. Univariate models:
- PD-L1 expression (stranded samples only)

2. Multivariate models:
- CD274 (PD-L1) expression, controlling for RNA library
- `quanTIseq immune cell fractions + extent of tumor resection + LGG broad histology`, across entire PBTA cohort
- `quanTIseq immune cell fractions + CD274 expression + extent of tumor resection + LGG broad histology` across entire PBTA cohort
- `quanTIseq immune cell fractions + CD274 expression + extent of tumor resection` across suitable `cancer_groups` (DECEASED N>=3, thus no LGG)


#### Output

Within the `results` folder are folders containing survival model fit results for each analysis described above.
```
results/immune/*.RDS (`survival-analysis_immune.Rmd`)
results/subtypes/*.RDS (`survival-analysis_subtypes.Rmd`)
results/template/*.tsv (`survival-analysis_template.Rmd`)
results/tp53_telomerase/*.RDS (`survival-analysis_tp53_telomerase.Rmd`)
```

The `plots` folder contains Kaplan-meier curve for the log rank analysis listed in parentheses below.
```
plots/KM_hgg_subtypes.pdf (`survival-analysis_subtypes.Rmd`)
plots/survival_curve_gender.pdf (`survival-analysis_template.Rmd`)
```
