# Batch Effect Analysis

Batch effects are commonly observed in RNA expression data. Such effects are likely less pronounced in RNA-Seq data than microarray data. However, batch effects may bias conclusions of studies that do not account for these effects. We wish to evaluate the level to which batch effects are present in the RNA-Seq data from this study and create alternative versions of the data that have been adjusted for batch effects.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Methods and Output](#methods-and-output)
 

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

This analysis can be run via the command line from the top directory of the
repository as follows:

```
bash run_batch_effects.sh
```

## Methods and Output

We tested for batch effects in three different areas:
1.	Batch = Sequencing Center
2.	Batch = Cohort
3.	Batch = polyA vs ribodepletion

Within each of these analyses I evaluated RNA seq data that had been prepared using either RSEM normalization or Kallisto (tpm) values.
The null hypothesis is that there are no batch effects in the data


|          |               | Batch                                          |                                          |                                    |
|----------|---------------|-------------------------------------------------------|------------------------------------------|------------------------------------|
|          |               | **Sequence center**                           | **Cohort**                               | **Preparation method**                 |
| RSEM     | polyA         | **reject**                                      | **reject**                                   | **reject**                             |
|          |               | before                                        | before                                   | before                             |
|          |               | - [rsem_polyA_sequence.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234696/rsem_polyA_sequence.pdf)                   | - [rsem_polya_cohort_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234692/rsem_polya_cohort_pca.pdf)            | - [rsem_method.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234686/rsem_method.pdf)               |
|          |               | - [rsem_polyA_sequence_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234687/rsem_polyA_sequence_pca.pdf)               | - [rsem_polyA_cohort.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234697/rsem_polyA_cohort.pdf)                | - [rsem_method_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234691/rsem_method_pca.pdf)            |
|          |               | after                                         | after                                    | after                              |
|          |               | - [rsem_polya_sequence_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325812/rsem_polya_sequence_combat.pdf)            | - [rsem_polya_cohort_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325813/rsem_polya_cohort_combat_pca.pdf)     | - [rsem_method_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325817/rsem_method_combat.pdf)        |
|          |               | - [rsem_polya_sequence_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325811/rsem_polya_sequence_combat_pca.pdf)        | - [rsem_polya_cohort_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325814/rsem_polya_cohort_combat.pdf)         | - [rsem_method_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325815/rsem_method_combat_pca.pdf)    |
|          |               |                                               |                                          |                                    |
|          | ribodepletion | **reject**                                        | **fail to reject**                           |                                    |
|          |               | before                                        | before                                   |                                    |
|          |               | - [rsem_stranded_sequence.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234683/rsem_stranded_sequence.pdf)               | - [rsem_stranded_cohort.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234695/rsem_stranded_cohort.pdf)             |                                    |
|          |               | - [rsem_stranded_sequence_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234681/rsem_stranded_sequence_pca.pdf)           | - [rsem_stranded_cohort_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234693/rsem_stranded_cohort_pca.pdf)         |                                    |
|          |               | after                                         | after                                    |                                    |
|          |               | - [rsem_stranded_sequence_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325809/rsem_stranded_sequence_combat_pca.pdf)     |                                    |                                    |
|          |               | - [rsem_stranded_sequence_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325810/rsem_stranded_sequence_combat.pdf)         |                                    |                                    |
|          |               |                                               |                                          |                                    |
| Kallisto | polyA         | **reject**                                        | **reject**                                   | **reject**                             |
|          |               | before                                        | before                                   | before                             |
|          |               | - [kallisto_polyA_sequence.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234682/kallisto_polyA_sequence.pdf)               | - [kallisto_polyA_cohort.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234684/kallisto_polyA_cohort.pdf)            | - [kallisto_method_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234694/kallisto_method_pca.pdf)        |
|          |               | - [kallisto_polya_sequence_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234680/kallisto_polya_sequence_pca.pdf)           | - [kallisto_polyA_cohort_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234685/kallisto_polyA_cohort_pca.pdf)        | - [kallisto_method.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234688/kallisto_method.pdf)            |
|          |               | after                                         | after                                    | after                              |
|          |               | - [kallisto_polya_sequence_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325819/kallisto_polya_sequence_combat_pca.pdf)    | - [kallisto_polya_cohort_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4328136/kallisto_polya_cohort_combat.pdf)     | - [kallisto_method_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325822/kallisto_method_combat.pdf)    |
|          |               | - [kallisto_polya_sequence_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325820/kallisto_polya_sequence_combat.pdf)        | - [kallisto_polya_cohort_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4328137/kallisto_polya_cohort_combat_pca.pdf) | - [kallisto_method_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325821/kallisto_method_combat_pca.pdf) |
|          | ribodepletion | **reject**                                        | **fail to reject**                           |                                    |
|          |               | before                                        | before                                   |                                    |
|          |               | - [kallisto_stranded_sequence_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234699/kallisto_stranded_sequence_pca.pdf)        | - [kallisto_stranded_cohort_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234689/kallisto_stranded_cohort_pca.pdf )    |                                    |
|          |               | - [kallisto_stranded_sequence.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234690/Kallisto_stranded_sequence.pdf)            | - [kallisto_stranded_cohort.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4234698/kallisto_stranded_cohort.pdf)        |                                    |
|          |               | after                                         | after                                    |                                    |
|          |               | - [kallisto_stranded_sequence_combat.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325818/kallisto_stranded_sequence_combat.pdf)      |                                    |                                    |
|          |               | - [kallisto_stranded_sequence_combat_pca.pdf](https://github.com/AlexsLemonade/OpenPBTA-analysis/files/4325816/kallisto_stranded_sequence_combat_pca.pdf)      |                                   |                                    |
|          |               |                                               |                                          |                                    | 

We propose to use the BatchQC tool to evaluate the data for batch effects. BatchQC provides various visualization and quantitative tools for evaluating batch effects. We will use these and prepare a summary based on our findings. ComBat is a widely used method that uses empirical Bayes methods for correcting batch effects. 

We used three covariates as batch variables in three separate experiments. 1) The scientific cohort 2) the sequencing center and 3) the library preparation method (polyA vs ribodepletion).



