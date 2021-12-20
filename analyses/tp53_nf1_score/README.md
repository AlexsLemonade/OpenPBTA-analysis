## Apply classifiers trained on TCGA RNA-seq data

**Module author:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6)) and Jaclyn Taroni (@jaclyn-taroni); code adapted from Gregory Way ([@gwaygenomics](https://github.com/gwaygenomics))

**Modified:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

This module is adapted from: [`marislab/pdx-classification`](https://github.com/marislab/pdx-classification).
Now published in [Rokita et al. _Cell Reports._ 2019.](https://doi.org/10.1016/j.celrep.2019.09.071)

In brief, _TP53_ inactivation, _NF1_ inactivation, and Ras activation classifiers are applied to the stranded and polya OpenPBTA RNA-seq data.
The classifiers were trained on TCGA PanCan data ([Way et al. _Cell Reports._ 2018](https://doi.org/10.1016/j.celrep.2018.03.046), [Knijnenburg et al. _Cell Reports._ 2018.](https://doi.org/10.1016/j.celrep.2018.03.076)).
See [`01-apply-classifier.py`](01-apply-classifier.py) for more information about the procedure.
To evaluate the classifier scores, we use [`06-evaluate-classifier.py`](06-evaluate-classifier.py) and input SNV data to identify true TP53/NF1 loss samples and compare scores of shuffled data to true calls and plot ROC curves. 

#### Running the analysis

The analysis can be run with the following (assuming you are in the root repository of the project):

```
bash analyses/tp53_nf1_score/run_classifier.sh
```

### Inputs from data download

* `pbta-snv-consensus-mutation.maf.tsv.gz`: from `snv-callers` module which gathers calls that are present in all 3 callers (strelka2,mutect2 and lancet) 
* `pbta-snv-scavenged-hotspots.maf.tsv.gz`: from `hotspot-detection` module to gather calls that overlap MSKCC hotspots found in any caller (except if only vardict calls the site as variant, we remove these calls since we have a lot of calls unique to vardict which we consider as false positive as discussed [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#snv-caller-comparison-analysis))
* pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds : stranded FPKM values per gene in biospecimen_id
* pbta-gene-expression-rsem-fpkm-collapsed.polya.rds : polya FPKM values per gene in biospecimen_id
* consensus_seg_with_status.tsv : created by analyses/focal-cn-file-preparation/02-add-ploidy-consensus.Rmd
* pbta-sv-manta.tsv.gz : Structural Variants called by manta


### Order of analysis

* `00-tp53-nf1-alterations.R` produces `TP53_NF1_snv_alteration.tsv`, which contains information about the presence or absence of coding SNVs in _TP53_ and _NF1_ for the purpose of evaluating the classifier results.
For evaluation purposes, a coding SNV 1) is found in a CDS region and 2) is not a silent mutation or in an intron as indicated by the `Variant_Classification` column of the consensus mutation file.
_NF1_ positive examples are additionally filtered to remove missense mutations, as these are not annotated with OncoKB ([#381 (comment)](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/381#issuecomment-570748578)).

* `01-apply-classifier.py` produces  `results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores.tsv`  and `results/pbta-gene-expression-rsem-fpkm-collapsed.polya_classifier_scores.tsv`, which contains all 3 classifier scores for the stranded data and for shuffled stranded (e.g., random) data.

* `02-qc-rna_expression_score.Rmd` here expression of TP53 gene was compared to TP53 classifier score. We didn't find a strong correlation between TP53 expression and TP53 inactivation score, thus, expression and classifier score together cannot predict function.

* `03-tp53-cnv-loss-domain.Rmd` here copy_number of regions overlapping TP53 functional domains were compared to TP53 classifier score to visualize the distribution of scores. We also save the TP53 loss calls which will be input for downstream altered status script.

* `04-tp53-sv-loss.Rmd` here structural variant breakpoints overlapping TP53 are investigated for CNV loss or low expression to gather high confidence TP53 loss via Structural Variants. 

* `05-tp53-altered-annotation.Rmd` here we take a deeper look into tp53 altered status with respect to number of SNVs/CNVs suggesting bi-allelic mutations or with respect to cancer_predisposition, molecular subtypes and tp53 classifier scores.


Columns | Description
-- | --
sample_id	| 7316-XXXX id used to match DNA and RNA Kids_First_Biospecimen_IDs
Kids_First_Biospecimen_ID_DNA	| Associated DNA Kids_First_Biospecimen_IDs
Kids_First_Biospecimen_ID_RNA	| Associated RNA Kids_First_Biospecimen_IDs
cancer_predispositions	| Germline cancer predisposition status
tp53_score	| TP53 loss classifier score
SNV_indel_counts	| Number of deleterious SNVs found in DNA sample
CNV_loss_counts	| Number of CNV losses found in DNA sample
HGVSp_Short	| Short format of protein level change used as SNV evidence 
CNV_loss_evidence | copy_number of CNV overlapping functional domains of TP53 used as evidence 	
hotspot	| Any 1 SNV shown in HGVSp_Short overlaps MSKCC cancer hotspot [database](https://www.cancerhotspots.org/#/home)
activating	| Any 1 SNV shown in HGVSp_Short overlaps TP53_ activating mutations R273C and R248W. [Reference](https://pubmed.ncbi.nlm.nih.gov/17417627/) and [reference](https://pubmed.ncbi.nlm.nih.gov/24677579/). 
tp53_altered | Combined evidence, cancer predisposition and score based tp53 status

**TP53 expression profile for activating vs loss TP53 status:**

```
plots
├── tp53_expression_by_altered_status_polya.png
├── tp53_expression_by_altered_status_stranded.png
```

**Plots of TP53 scores vs TP53 altered status:**

```
plots
├── tp53_scores_by_altered_status.png
```

**Plots of TP53 scores vs TP53 altered status per cancer predispositions:**

```
plots
├── tp53_scores_vs_tp53_altered_status_Li-Fraumeni\ syndrome.png
├── tp53_scores_vs_tp53_altered_status_NF-1,Other\ inherited\ conditions\ NOS.png
├── tp53_scores_vs_tp53_altered_status_NF-1.png
├── tp53_scores_vs_tp53_altered_status_NF-2,Schwannomatosis.png
├── tp53_scores_vs_tp53_altered_status_NF-2.png
├── tp53_scores_vs_tp53_altered_status_None\ documented.png
├── tp53_scores_vs_tp53_altered_status_Other\ inherited\ conditions\ NOS,Schwannomatosis.png
├── tp53_scores_vs_tp53_altered_status_Other\ inherited\ conditions\ NOS.png
├── tp53_scores_vs_tp53_altered_status_Schwannomatosis.png
└── tp53_scores_vs_tp53_altered_status_Tuberous\ Sclerosis.png
```

* `06-evaluate-classifier.py` evaluates classifier score with TP53 alterations (non-synonymous SNV and all status == "loss" in consensus CNV file from 00-tp53-nf1-alterations.R) 

**ROC threshold results for shuffled vs non-shuffled stranded and polya classifier output:**

```
results
├── polya_TP53_roc_threshold_results.tsv
├── polya_TP53_roc_threshold_results_shuffled.tsv
├── stranded_TP53_roc_threshold_results.tsv
├── stranded_TP53_roc_threshold_results_shuffled.tsv
```

Because some of the classifier genes are not present in the OpenPBTA dataset, the scores should be interpreted as continuous values representing relative gene alterations and not as probabilities.

* `07-plot-roc.R` using the output obtained from `06-evaluate-classifier.py`, it creates ROC curves for TP53 classifier scores for stranded and poly-A RNA-seq data: 

```
plots
├── polya_TP53_roc.png
├── stranded_TP53_roc.png
```

* `08-compare-molecularsubtypes-tp53scores.R` creates violin plots of TP53 scores across all molecular subtypes per broad histology.

**Plots of TP53 scores vs Molecular subtypes per broad histology:**

```
plots
├── tp53_scores_vs_molecular_subtype_Diffuse_astrocytic_and_oligodendroglial_tumor.png
├── tp53_scores_vs_molecular_subtype_Embryonal_tumor.png
├── tp53_scores_vs_molecular_subtype_Ependymal_tumor.png
├── tp53_scores_vs_molecular_subtype_Low-grade_astrocytic_tumor.png
```

**Associated p-values:**

```
results
├── tp53_scores_vs_molecular_subtype_Diffuse_astrocytic_and_oligodendroglial_tumor.tsv
├── tp53_scores_vs_molecular_subtype_Embryonal_tumor.tsv
├── tp53_scores_vs_molecular_subtype_Ependymal_tumor.tsv
└── tp53_scores_vs_molecular_subtype_Low-grade_astrocytic_tumor.tsv
```

* `09-compare-histologies.R` creates boxplots comparing broad histology and cancer group.

**Boxplots:**

```
plots
├── tp53_broad_histology.pdf
├── tp53_cancer_group.pdf
```
