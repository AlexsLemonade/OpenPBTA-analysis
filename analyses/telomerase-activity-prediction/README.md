Telomerase Activity Prediction

Module Authors: Nighat Noureen and Siyuan Zheng. 

Goals include:

1. Predict telomerase activities using gene expression data.

2. Analyze the variations in telomerase activities across polyA and stranded data sets.

3. Analyze the correlation of telomerase activities with TERT and TERC genes expression across polyA and stranded data sets.

4. Analyze the distribution of telomerase activities across different histologies (used short and broad histologies from the clinical data)

5. Generate plots for the analysis modules.




Scripts and Results/Plots:

1. RUN-telomerase-activity-prediction.sh is used to generate the telomerase activities from gene expression data sets from different platforms.

2. Results directory contains Telomerase Scores for PolyA counts data (file = "TelomeraseScores_PTBAPolya_counts.txt"),PolyA FPKM data (file = "TelomeraseScores_PTBAPolya_FPKM.txt"), Stranded Counts data(file = "TelomeraseScores_PTBAStranded_counts.txt") and Stranded FPKM data (file = "TelomeraseScores_PTBAStranded_FPKM.txt")

3. Analysis-telomerase-activity.sh is used to perform comparative analysis and the results are saved in plots directory.

4. "PTBA_GE_Score_AllScatter.pdf" plot represents the telomerase activity scores correlations across different data sets (i.e Stranded counts versus Stranded FPKM and PolyA counts versus PolyA FPKM)

5. "PTBA_GE_TM_ScatterComp.pdf" plot correlates TERT and TERC gene expressions with telomerase activities for different data sets (i.e Stranded and PolyA)

6. "PBTA_StrandedHistology.pdf" plot shows the distribution of telomerase activities across different brain tumor histologies.