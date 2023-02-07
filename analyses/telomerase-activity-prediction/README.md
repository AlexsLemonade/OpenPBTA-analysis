Telomerase Activity Prediction

Module Authors: Nighat Noureen and Siyuan Zheng.

Goals include:

1. Predict telomerase activities using gene expression data.

2. Analyze the variations in telomerase activities across polyA and stranded data sets.

3. Analyze the correlation of telomerase activities with TERT and TERC genes expression across polyA and stranded data sets.

4. Analyze the distribution of telomerase activities across different histologies (used short and broad histologies from the clinical data).

5. Compare the telomerase activities across different molecular subtypes of various PBTA histologies (used short histologies from the clinical data).

6. Generate plots for the analysis modules.

7. Explore how results for stranded libraries may be influenced if only tumors passing a given tumor purity threshold are considered.


To run the complete module run as follows (assuming you are in this directory):
```sh
bash RUN-telomerase-activity-prediction.sh
```


To re-run generate analysis results that are presented in the manuscript, run as follows (again, assuming you are in this directory):

```sh
 OPENPBTA_FOR_FIGURES=1 bash RUN-telomerase-activity-prediction.sh
```



Contents of directory:

1. [`RUN-telomerase-activity-prediction.sh`](./RUN-telomerase-activity-prediction.sh) is used to generate the telomerase activities from gene expression data sets from different platforms, and calls all R scripts in order.


2. The directory `results/` contains:
	+ Telomerase Scores for PolyA counts data in [TelomeraseScores_PTBAPolya_counts.txt](./results/TelomeraseScores_PTBAPolya_counts.txt)
	+ Telomerase Scores for PolyA FPKM data in [TelomeraseScores_PTBAPolya_FPKM.txt](./results/TelomeraseScores_PTBAPolya_FPKM.txt)
	+ Telomerase Scores for Stranded Counts data in [TelomeraseScores_PTBAStranded_counts.txt](./results/TelomeraseScores_PTBAStranded_counts.txt)
	+ Telomerase Scores for Stranded FPKM data in [TelomeraseScores_PTBAStranded_FPKM.txt](./results/TelomeraseScores_PTBAStranded_FPKM.txt)
	+ Telomerase Scores for Stranded FPKM data filtered to a cancer-group-specific tumor purity threshold in in [TelomeraseScores_PTBAStranded_FPKM_thresholded.txt](./results/TelomeraseScores_PTBAStranded_FPKM_thresholded.txt)
	+ `EXTENDScores_{broad_histology}.txt`: compares telomerase activities across different molecular subtypes of various PBTA short histologies and contains p-values as well as adjusted p-values.

3. The directory `plots/` contains:
	+ [PTBA_GE_Score_AllScatter.pdf](./plots/PTBA_GE_Score_AllScatter.pdf): represents the telomerase activity scores correlations across different data sets (i.e Stranded counts versus Stranded FPKM and PolyA counts versus PolyA FPKM)
	+ [PTBA_GE_TM_ScatterComp.pdf](./plots/PTBA_GE_TM_ScatterComp.pdf): correlates TERT and TERC gene expressions with telomerase activities for different data sets (i.e Stranded and PolyA)
	+ [PBTA_StrandedHistology.pdf](./plots/PBTA_StrandedHistology.pdf): shows the distribution of telomerase activities across different brain tumor histologies
	+ `EXTENDScores_{broad_histology}.png`: compares telomerase activities across different molecular subtypes of various PBTA short histologies.
	+ [`TERTp_mutations.pdf`](./plots/TERTp_mutations.pdf): EXTEND score distribition is shown for samples with and without specific TERT promoter mutations.
	+ [`TERTp_mutations_TERC_TERT_expression.pdf`](./plots/TERTp_mutations_TERC_TERT_expression.pdf): Relationship between TERT and TERC stranded expression and EXTEND scores, highlighting samples with TERTp mutations
