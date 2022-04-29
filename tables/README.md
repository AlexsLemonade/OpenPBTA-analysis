## Tables Generation for paper

This module creates tables and supplementary tables for the manuscript 

Table 1: Molecular subtypes determined for this project
- This table shows the type and number of different molecular subtypes profiled in the project.

Table 2: Patients with hypermutant tumors
- Listed are patients with at least one hypermutant or ultra-hypermutant tumor (inclusive of derived cell lines and all phases of therapy).
Coding region TMB, phase of therapy, therapeutic interventions, cancer predispositions, and molecular subtypes are included.

Table S1: V21 histologies table
- This table displays the histologies file from v21 release of OpenPBTA in an excel format.
    - Sheet 1: README
    - Sheet 2: V21 histologies file
    - Sheet 3: CNS region definitions

Table S2: DNA results table
- This excel shows key results from analyzing DNA sequencing. 
    - Sheet 1: this table shows tumor mutation burden (TMB) for all samples - `TMB_coding` is calculated based on mutations in coding regions and `TMB_all` is calculated by all mutations a sample has
        - Input files `pbta-snv-consensus-mutation-tmb-all.tsv` and `pbta-snv-consensus-mutation-tmb-coding.tsv` are in data release.
    - Sheet 2: this table shows signature weights for CNS signature identified in [Degasperi et al, 2020](https://doi.org/10.1038/s43018-020-0027-5) using [`deconstructSigs`](https://doi.org/10.1186/s13059-016-0893-4) per sample
        - Input file `deconstructsigs_exposures_merged.tsv` is generated from module [mutational-signatures](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)
    - Sheet 3: this table summarizes the number of chromothripsis events per sample 
        - Input file `chromothripsis_summary_per_sample.txt` is generated from module [chromothripsis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis)
    
Table S3: RNA results table
- This excel shows key results from analyzing RNA sequencing. 
    - Sheet 1: this table shows TP53 scores, cancer prediposition, SNV indel counts, CNV loss counts, SV counts, fusion counts and details about these events for each sample 
        - Input file `tp53_altered_status.tsv` is from module [tp53_nf1_scores](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tp53_nf1_score)
    - Sheet 2: this tables shows the EXTEND scores for all samples - `NormEXTENDScores_counts` is EXTEND scores calculated based on expression counts matrix and `NormEXTENDScores_fpkm` is EXTEND scores calculated based on expression FPKM matrix
        - Input files are from module [telomerase-activity-prediction](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction)
    - Sheet 3: this tables shows the fractions of for various immune cell types in the tumor microenvironment (TME) for each RNA sample as calculated with [`quanTIseq`](https://doi.org/10.1186/s13073-019-0638-6).
        - Input files are from module [immune-deconv](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv)
    
Table S4: Key Resources
- Table of all software and their respective versions used for the OpenPBTA project. Of note, this table contains all software in the OpenPBTA docker image utilized within the repository, but not all software was used for the final manuscript.
    - Sheet 1: r_packages
    - Sheet 2: python_libraries
    - Sheet 3: other_command_line_tools
    - Sheet 4: workflow_repository_tools

    
