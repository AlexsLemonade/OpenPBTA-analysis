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
library(forcats) ### for fct_reorder

############################################################### Combining Histology with EXTEND Scores (Figure 3) ###############################################################################



PTBA_GE_Standard_TMScores = read.table('results/TelomeraseScores_PTBAStranded_FPKM.txt',sep='\t',head=T)

PTBA_Histology = read.table('../../data/pbta-histologies.tsv',sep='\t',head=T)
colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'
PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')

Stranded_Histology = PTBA_GE_Standard_Histology
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology == "Other"),]
Frequency = data.frame(table(Stranded_Histology$short_histology))
colnames(Frequency)=c('Variables','Freq')
Frequency = Frequency[which(Frequency$Freq == 1),]
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology %in% Frequency$Variables),]


pdf('PBTA_StrandedHistology.pdf',onefile=FALSE)

P1 = ggplot(Stranded_Histology, aes(x=fct_reorder(short_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="cyan3")+theme_classic()+theme(axis.text.x=element_text(angle=50,size=7,vjust=1,hjust=1),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))

P2 = ggplot(Stranded_Histology, aes(x=fct_reorder(broad_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size= 0.1,notch=FALSE,outlier.size = 0,outlier.shape=NA,fill="cyan3")+theme_classic()+theme(axis.text.x=element_text(angle=50,size=7,vjust=1,hjust=1),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))


grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =4, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",6),font.ytickslab =c("black",6),font.x = 7,font.y=7,font.legend=6,xlab="Tumor Histology(short)",ylab="EXTEND Scores"),vp = define_region(row = 1:2, col = 1:2))
print(ggpar(P2,font.xtickslab =c("black",6),font.ytickslab =c("black",6),font.x = 7,font.y=7,font.legend=6,xlab="Tumor Histology(broad)",ylab="EXTEND Scores"),vp = define_region(row = 3:4, col = 1:2))


dev.off()


