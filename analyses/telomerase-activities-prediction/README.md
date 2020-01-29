Telomerase Activity Prediction

Module Authors: Nighat Noureen and Siyuan Zheng. 

Goals include:

1. Predict telomerase activities using gene expression data.

2. Analyze the variations in telomerase activities across polyA and stranded data sets.

3. Analyze the correlation of telomerase activities with TERT and TERC genes expression across polyA and stranded data sets.

4. Analyze the distribution of telomerase activities across different histologies (used short and broad histologies from the clinical data)

5. Generate plots for the analysis modules.




Scripts and Results:

1. 01-RSEM-FPKM-Stranded.r calculates telomerase activities for RSEM-FPKM-Stranded data. Result file is saved in results/TelomeraseScores_PTBAStranded_FPKM.txt file(Data V13).

2. 02-RSEM-FPKM-PolyA.r calculates telomerase activities for RSEM-FPKM-PolyA data. Result file is saved in results/TelomeraseScores_PTBAPolyA_FPKM.txt file(Data V13).

3. 03-RSEM-Counts-Stranded.r  calculates telomerase activities for RSEM-Counts-Stranded data. Result file is saved in results/TelomeraseScores_Stranded_RsemCounts.txt file(Data V13).

4. 04-RSEM-Counts-PolyA.r  calculates telomerase activities for RSEM-Counts-PolyA data. Result file is saved in results/TelomeraseScores_PolyA_Rsemcounts.txt file(Data V13).

5. 05-Comparing-Counts-versus-FPKM.r  compares telomerase activities using counts and FPKM cases across polyA and stranded data sets. Result file is saved in plots/PTBA_GE_Score_AllScatter.pdf file(Data V13).

6. 06-Comparing-TERTexp-TERCexp-EXTENDScores.r compares TERT and TERC genes expression using counts and FPKM data with telomerase activities (EXTENDScores) across PolyA and Stranded data sets. Result file is saved in plots/PTBA_GE_TM_ScatterComp.pdf file(Data V13).

7. 07-Comparing-Histology-versus-EXTENDScores.r distributes the telomerase activities (EXTENDScores) across different histologies (using broad and short histologies). Result file is saved in plots/PBTA_StrandedHistology.pdf file(Data V13).
