## Tables Generation for paper

This module creates tables and supplementary tables for the manuscript 

Table 1: Molecular subtypes determined for this project
- This table shows the type and number of different molecular subtypes profiled in the project.

Table S1: V21 histologies table
- This table displays the histologies file from v21 release of OpenPBTA in an excel format.

Table S2: DNA results table
- This excel file aims to include key results from analyzing DNA sequencing results 
    - Sheet 1: this table shows tumor mutation burden (TMB) for all samples - `TMB_coding` is calculated based on mutations in coding regions and `TMB_all` is calculated by all mutations a sample has
        - Input files `pbta-snv-mutation-tmb-all.tsv` and `pbta-snv-mutation-tmb-coding.tsv` are in data release.
    - Sheet 2: this table shows signature weights for CNS signature identified in [Degasperi et al, 2020](https://doi.org/10.1038/s43018-020-0027-5) using [`sigfit`](https://github.com/kgori/sigfit) per sample
        - Input file `fitted_exposures_signal-cns_sigfit.rds` is generated from module [mutational-signatures](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)
    - Sheet 3: this table summarizes the number of chromothripsis events per sample 
        - Input file `chromothripsis_summary_per_sample.txt` is generated from module [chromothripsis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis)
    
Table S3: RNA results table
- This excel file aims to include key results from analyzing RNA sequencing results 
    - Sheet 1: this table shows TP53 scores, cancer prediposition, SNV indel counts, CNV loss counts, SV counts, fusion counts and details about these events for each sample 
        - Input file `tp53_altered_status.tsv` is from module [tp53_nf1_scores](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tp53_nf1_score)
    - Sheet 2: this tables shows the EXTEND scores for all samples - `NormEXTENDScores_counts` is EXTEND scores calculated based on expression counts matrix and `NormEXTENDScores_fpkm` is EXTEND scores calculated based on expression FPKM matrix
        - Input files are from module [telomerase-activity-prediction](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction)
    - Sheet 3: this tables shows quantifies scores for various immune cell types in the tumor microenvironment (TME) for each RNA sample -
        - Input files are from module [immune-deconv](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv)
    
    
    
