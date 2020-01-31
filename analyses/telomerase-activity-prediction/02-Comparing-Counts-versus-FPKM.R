library('pracma')
library('ggplot2')
library('gplots')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(gridBase)
library(gridGraphics)




################################################ Comparing Counts versus FPKM (Figure 1) ############################################################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_counts.txt")   ### Variable representing Telomerase activity-prediction for Polya_counts
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_counts.txt")  ### Variable representing Telomerase activity-prediction for Stranded_counts
				 
C = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_FPKM.txt")      ### Variable representing Telomerase activity-prediction for PolyA_FPKM data
				 
D = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data

PTBA_GE_PolyA_TMScores1 = read.table(A,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores1)[colnames(PTBA_GE_PolyA_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_PolyACounts'

PTBA_GE_Standard_TMScores1 = read.table(B,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores1)[colnames(PTBA_GE_Standard_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_StrandedCounts'

PTBA_GE_PolyA_TMScores2 = read.table(C,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores2)[colnames(PTBA_GE_PolyA_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_PolyA_FPKM'

PTBA_GE_Standard_TMScores2 = read.table(D,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores2)[colnames(PTBA_GE_Standard_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_Stranded_FPKM'

PTBA_GE_PolyA_TMScores = merge(PTBA_GE_PolyA_TMScores1,PTBA_GE_PolyA_TMScores2,by='SampleID')
PTBA_GE_Standard_TMScores = merge(PTBA_GE_Standard_TMScores1,PTBA_GE_Standard_TMScores2,by='SampleID')


pdf('plots/PTBA_GE_Score_AllScatter.pdf')


P1 = ggscatter(PTBA_GE_PolyA_TMScores, x = "NormEXTENDScores_PolyACounts", y = "NormEXTENDScores_PolyA_FPKM",color = "red",size = 2,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)

		   
P2 = ggscatter(PTBA_GE_Standard_TMScores, x = "NormEXTENDScores_StrandedCounts", y = "NormEXTENDScores_Stranded_FPKM",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)
		   
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND_PolyA_Counts",ylab="EXTEND_PolyA_FPKM"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND_Stranded_Counts",ylab="EXTEND_Stranded_FPKM"),vp = define_region(row = 1, col = 2))

		   
dev.off()

