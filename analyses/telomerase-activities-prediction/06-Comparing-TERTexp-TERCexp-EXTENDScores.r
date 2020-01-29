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




###############################################  Comparing TERT and TERC expression with EXTEND Scores (Figure 2)  ################################################################################################



PTBA_GE_PolyA = readRDS('../../data/pbta-gene-expression-rsem-fpkm.polya.rds')
Genes = substr(PTBA_GE_PolyA$gene_id, stringr::str_locate(PTBA_GE_PolyA$gene_id,"_")+1,stringr::str_length(PTBA_GE_PolyA$gene_id))
PTBA_GE_PolyA = data.frame(Genes,PTBA_GE_PolyA)
PTBA_GE_PolyA = PTBA_GE_PolyA[,-2]
TelComp = c('TERT','TERC')
PTBA_GE_PolyA_GE = PTBA_GE_PolyA[which(PTBA_GE_PolyA$Genes %in% TelComp),]

PTBA_GE_Standard = readRDS('../../data/pbta-gene-expression-rsem-fpkm.stranded.rds')
Genes = substr(PTBA_GE_Standard$gene_id, stringr::str_locate(PTBA_GE_Standard$gene_id,"_")+1,stringr::str_length(PTBA_GE_Standard$gene_id))
PTBA_GE_Standard = data.frame(Genes,PTBA_GE_Standard)
PTBA_GE_Standard = PTBA_GE_Standard[,-2]
TelComp = c('TERT','TERC')
PTBA_GE_Standard_GE = PTBA_GE_Standard[which(PTBA_GE_Standard$Genes %in% TelComp),]


PTBA_GE_PolyA_TMScores2 = read.table('results/TelomeraseScores_PTBAPolyA_FPKM.txt',sep='\t',head=T)
PTBA_GE_Standard_TMScores2 = read.table('results/TelomeraseScores_PTBAStranded_FPKM.txt',sep='\t',head=T)


rownames(PTBA_GE_PolyA_GE)= PTBA_GE_PolyA_GE$Genes
PTBA_GE_PolyA_GE = PTBA_GE_PolyA_GE[,-1]
PTBA_GE_PolyA_GE = t(PTBA_GE_PolyA_GE)
SampleID = rownames(PTBA_GE_PolyA_GE)
PTBA_GE_PolyA_GE = data.frame(SampleID,PTBA_GE_PolyA_GE)
row.names(PTBA_GE_PolyA_GE)=NULL


rownames(PTBA_GE_Standard_GE)= PTBA_GE_Standard_GE$Genes
PTBA_GE_Standard_GE = PTBA_GE_Standard_GE[,-1]
PTBA_GE_Standard_GE = t(PTBA_GE_Standard_GE)
SampleID = rownames(PTBA_GE_Standard_GE)
PTBA_GE_Standard_GE = data.frame(SampleID,PTBA_GE_Standard_GE)
row.names(PTBA_GE_Standard_GE)=NULL



PTBA_GE_TM_PolyA = merge(PTBA_GE_PolyA_TMScores2,PTBA_GE_PolyA_GE,by='SampleID')
PTBA_GE_TM_Standard = merge(PTBA_GE_Standard_TMScores2,PTBA_GE_Standard_GE,by='SampleID')



pdf('PTBA_GE_TM_ScatterComp.pdf')


P1 = ggscatter(PTBA_GE_TM_PolyA, x = "NormEXTENDScores", y = "TERT",color = "red",size = 2,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)

		   
P2 = ggscatter(PTBA_GE_TM_PolyA, x = "NormEXTENDScores", y = "TERC",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)
		   
P3 = ggscatter(PTBA_GE_TM_Standard, x = "NormEXTENDScores", y = "TERT",color = "red",size = 1.3,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)
		   	   

P4 = ggscatter(PTBA_GE_TM_Standard, x = "NormEXTENDScores", y = "TERC",color = "red",size = 1.3,
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


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores PolyA_FPKM",ylab="TERT Expression"),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores PolyA_FPKM",ylab="TERC Expression"),vp = define_region(row = 1, col = 2))
print(ggpar(P3,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores Stranded_FPKM",ylab="TERT Expression"),vp = define_region(row = 2, col = 1))
print(ggpar(P4,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,xlab="EXTEND Scores Stranded_FPKM",ylab="TERC Expression"),vp = define_region(row = 2, col = 2))
		   
dev.off()

