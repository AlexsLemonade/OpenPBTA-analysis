library('ggpubr')
stringsAsFactors=FALSE
suppressPackageStartupMessages({
  library(stringr)
  library(gridBase)
  library(gridGraphics)
  library(optparse)
})

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


# set up the command line options
option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Plot output (.pdf)")
)


# parse the command line arguments - this becomes a list
# where each element is named by the long flag option
opt <- parse_args(OptionParser(option_list = option_list))

PBTA_GEComparisonPlot <- opt$output

				 
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


pdf(PBTA_GEComparisonPlot,height=2,width=4)


P1 = ggscatter(PTBA_GE_PolyA_TMScores, x = "NormEXTENDScores_PolyACounts", y = "NormEXTENDScores_PolyA_FPKM",color = "red",size = 0.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=1)

		   
P2 = ggscatter(PTBA_GE_Standard_TMScores, x = "NormEXTENDScores_StrandedCounts", y = "NormEXTENDScores_Stranded_FPKM",color = "red",size = 0.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=1)
		   
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =4,font.ytickslab =4,font.x = 4,font.y=4,xlab="EXTEND_PolyA_Counts",ylab="EXTEND_PolyA_FPKM"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =4,font.ytickslab =4,font.x = 4,font.y=4,xlab="EXTEND_Stranded_Counts",ylab="EXTEND_Stranded_FPKM"),vp = define_region(row = 1, col = 2))

		   
dev.off()

