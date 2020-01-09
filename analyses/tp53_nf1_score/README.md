## Apply classifiers trained on TCGA RNA-seq data

**Module author:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6))

This module is adapted from: [`marislab/pdx-classification`](https://github.com/marislab/pdx-classification).
Now published in [Rokita et al. _Cell Reports._ 2019.](https://doi.org/10.1016/j.celrep.2019.09.071)

In brief, _TP53_ inactivation, _NF1_ inactivation, and Ras activation classifiers are applied to the stranded and polya OpenPBTA RNA-seq data.
The classifiers were trained on TCGA PanCan data ([Way et al. _Cell Reports._ 2018](https://doi.org/10.1016/j.celrep.2018.03.046), [Knijnenburg et al. _Cell Reports._ 2018.](https://doi.org/10.1016/j.celrep.2018.03.076)).
See [`01-apply-classifier.py`](01-apply-classifier.py) for more information about the procedure.
To evaluate the classifier scores, we use [`02-evaluate-classifier.py`](02-evaluate-classifier.py) and input SNV data to identify true TP53/NF1 loss samples and compare scores of shuffled data to true calls and plot ROC curves. 

**This module is in progress.** 
Copy number aberrations are not currently considered when evaluating the classifiers.

#### Running the analysis

The analysis can be run with the following (assuming you are in the root repository of the project):

```
bash analyses/tp53_nf1_score/run_classifier.sh
```

#### Output

`00-tp53-nf1-alterations.R` produces `TP53_NF1_snv_alteration.tsv`, which contains information about the presence or absence of coding SNVs in _TP53_ and _NF1_ for the purpose of evaluating the classifier results.
For evaluation purposes, a coding SNV 1) is found in a CDS region and 2) is not a silent mutation or in an intron as indicated by the `Variant_Classification` column of the consensus mutation file.
_NF1_ positive examples are additionally filtered to remove missense mutations, as these are not annotated with OncoKB ([#381 (comment)](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/381#issuecomment-570748578)).

`01-apply-classifier.py` produces  `results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores.tsv`  and `results/pbta-gene-expression-rsem-fpkm-collapsed.polya_classifier_scores.tsv`, which contains all 3 classifier scores for the stranded data and for shuffled stranded (e.g., random) data.

Because some of the classifier genes are not present in the OpenPBTA dataset, the scores should be interpreted as continuous values representing relative gene alterations and not as probabilities.

ROC curve for TP53 classifier scores for stranded RNAseq data
![stranded RNAseq TP53 classifier ROC](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/validation_step/analyses/tp53_nf1_score/results/stranded_TP53.png)

ROC curve for TP53 classifier scores for polya RNAseq data
![polya RNAseq TP53 classifier ROC](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/validation_step/analyses/tp53_nf1_score/results/polya_TP53.png)

ROC curve for NF1 classifier scores for stranded RNAseq data
![stranded RNAseq NF1 classifier ROC](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/validation_step/analyses/tp53_nf1_score/results/stranded_NF1.png)

ROC curve for NF1 classifier scores for polya RNASeq data
![polya RNASeq NF1 classifier ROC](https://github.com/kgaonkar6/OpenPBTA-analysis/blob/validation_step/analyses/tp53_nf1_score/results/polya_NF1.png)


