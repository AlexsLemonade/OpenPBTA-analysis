library('pracma')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library('devtools')
library(EXTEND)

###################################################   Running Analysis on RSEM-FPKM data formats  #####################################################################

PTBA_GE_PolyA = readRDS('../../data/pbta-gene-expression-rsem-fpkm.polya.rds')
MAX <- apply(PTBA_GE_PolyA[-1], FUN=max, MARGIN=1)
PTBA_GE_PolyA = data.frame(MAX,PTBA_GE_PolyA)
Genes = substr(PTBA_GE_PolyA$gene_id, stringr::str_locate(PTBA_GE_PolyA$gene_id,"_")+1,stringr::str_length(PTBA_GE_PolyA$gene_id))
PTBA_GE_PolyA = data.frame(Genes,PTBA_GE_PolyA)
PTBA_GE_PolyA <- PTBA_GE_PolyA[order(-PTBA_GE_PolyA[,"MAX"]),] #Sorting the Genes with MAX value (Asced to Desc)
PTBA_GE_PolyA = PTBA_GE_PolyA[!duplicated(PTBA_GE_PolyA$Genes),]
row.names(PTBA_GE_PolyA)= PTBA_GE_PolyA$Genes
PTBA_GE_PolyA = PTBA_GE_PolyA[,-c(1,2,3)]
data = PTBA_GE_PolyA
data = as.matrix(data)
RunEXTEND(data)
file.rename((pattern="TelomeraseScores.txt"), paste0("TelomeraseScores_PTBAPolyA_FPKM.txt"))



