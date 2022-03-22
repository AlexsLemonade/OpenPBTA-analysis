## Perform survival analysis for samples with 

**Author of code and documentation:** [@runjin326](https://github.com/runjin326)

In this module, we perform survival anlayses by Kaplan-Meier, log-rank or cox regression for samples in OpenPBTA, using different covariates. Detailed information about covariates and methods used can be found in module contents.


### Usage
```sh
bash run-survival.sh
```

### Module contents

`survival-analysis_tempalte.Rmd` showed the usage of function `survival_analysis` to perform Kaplan-Meier, log-rank and cox regression survival analyses. 
Particularly, it demonstrated runninng Kaplan-Meier and log-rank survival analysis on categorial variable, using `germline_sex_estimate` as an example.
It also demonstrated how to run cox regression survival analysis on continuous variable, using `tmb` (tumor mutation burden) as an example.

`survival-analysis_HGG_DMG.Rmd` runs Kaplan-Meier survival analysis on two groups of samples: HGG, H3 wildtype and DMG, H3 K28. 
Two visualizations were generated for the same analyis. 

`survival-analysis_histology.Rmd` run the following analysis:
1. Univariate analysis 
a) cox regression method
- TP53 classifier score (as a continuous variable)
- EXTEND score (as a continuous variable)
- HGG vs. non-HGG
b) log rank method
- HGG vs. non-HGG

2. Multivariate analysis - all using cox regression
a) overall comparison with 3 covariates
- TP53 classifier score, EXTEND score and HGG vs. non-HGG
b) sub-population comparison with 2 covariates
- HGG group: TP53 classifier score and EXTEND score
- non-HGG group: TP53 classifier score and EXTEND score

3. Plots
- Bivariate distribution of TP53 classifier score and EXTEND score in density plot stratified by HGAT status
