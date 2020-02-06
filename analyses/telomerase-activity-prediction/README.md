Telomerase Activity Prediction

Module Authors: Nighat Noureen and Siyuan Zheng. 

Goals include:

1. Predict telomerase activities using gene expression data.

2. Analyze the variations in telomerase activities across polyA and stranded data sets.

3. Analyze the correlation of telomerase activities with TERT and TERC genes expression across polyA and stranded data sets.

4. Analyze the distribution of telomerase activities across different histologies (used short and broad histologies from the clinical data)

5. Generate plots for the analysis modules.




Contents of directory:

1. [`RUN-telomerase-activity-prediction.sh`](./RUN-telomerase-activity-prediction.sh) is used to generate the telomerase activities from gene expression data sets from different platforms, and calls all R scripts in order.

2. The directory `results/` contains:
	+ Telomerase Scores for PolyA counts data in [TelomeraseScores_PTBAPolya_counts.txt](./results/TelomeraseScores_PTBAPolya_counts.txt)
	+ Telomerase Scores for PolyA FPKM data in [TelomeraseScores_PTBAPolya_FPKM.txt](./results/TelomeraseScores_PTBAPolya_FPKM.txt)
	+ Telomerase Scores for Stranded Counts data in [TelomeraseScores_PTBAStranded_counts.txt](./results/TelomeraseScores_PTBAStranded_counts.txt)
	+ Telomerase Scores for Stranded FPKM data in [TelomeraseScores_PTBAStranded_FPKM.txt](./results/TelomeraseScores_PTBAStranded_FPKM.txt)

3. The directory `plots/` contains:
	+ [PTBA_GE_Score_AllScatter.pdf](./plots/PTBA_GE_Score_AllScatter.pdf): represents the telomerase activity scores correlations across different data sets (i.e Stranded counts versus Stranded FPKM and PolyA counts versus PolyA FPKM)
	+ [PTBA_GE_TM_ScatterComp.pdf](./plots/PTBA_GE_TM_ScatterComp.pdf): correlates TERT and TERC gene expressions with telomerase activities for different data sets (i.e Stranded and PolyA)
	+ [PBTA_StrandedHistology.pdf](./plots/PBTA_StrandedHistology.pdf): shows the distribution of telomerase activities across different brain tumor histologies.