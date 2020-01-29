library('pracma')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library('devtools')
install_github('NNoureen/EXTEND')
library(EXTEND)


###################################################   Running Analysis on RSEM-FPKM data formats  #####################################################################



PTBA_GE_Standard = readRDS('../../data/pbta-gene-expression-rsem-fpkm.stranded.rds')
MAX <- apply(PTBA_GE_Standard[-1], FUN=max, MARGIN=1)
PTBA_GE_Standard = data.frame(MAX,PTBA_GE_Standard)
Genes = substr(PTBA_GE_Standard$gene_id, stringr::str_locate(PTBA_GE_Standard$gene_id,"_")+1,stringr::str_length(PTBA_GE_Standard$gene_id))
PTBA_GE_Standard = data.frame(Genes,PTBA_GE_Standard)
PTBA_GE_Standard <- PTBA_GE_Standard[order(-PTBA_GE_Standard[,"MAX"]),] #Sorting the Genes with MAX value (Asced to Desc)
PTBA_GE_Standard = PTBA_GE_Standard[!duplicated(PTBA_GE_Standard$Genes),]
row.names(PTBA_GE_Standard)= PTBA_GE_Standard$Genes  ##
PTBA_GE_Standard = PTBA_GE_Standard[,-c(1,2,3)]
data = PTBA_GE_Standard
data = as.matrix(data)
RunEXTEND(data)#####EXTEND 
file.rename((pattern="TelomeraseScores.txt"), paste0("TelomeraseScores_PTBAStranded_FPKM.txt"))




