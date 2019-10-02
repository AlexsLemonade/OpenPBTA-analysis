##########################################
#Purpose: Code to run ssGSEA analysis
#Author: Pichai Raman
#Date: 9/23/2019
##########################################


#Call libraries
library("GSVA")
library("GSEABase");


#Read in data 
expData <- readRDS("../../data/pbta-gene-counts-rsem-expected_count.stranded.rds");

#Format RNA-Seq data
#Get to just one gene symbol per row
getGene <- function(x)
{
    strsplit(x, split="_")[[1]][2]
}
expData[,"max"] <- apply(expData[-1], FUN=max, MARGIN=1);
expData[,"Gene"] <- sapply(as.character(expData[,1]), FUN=getGene);
expData <- expData[order(-expData[,"max"]),]
expData <- expData[!duplicated(expData[,"Gene"]),]
rownames(expData) <- expData[,"Gene"]
expData <- expData[expData[,"max"]>0,]
expData <- dplyr::select(expData, dplyr::starts_with("BS_"))
expData <- expData[-1];
expData <- as.matrix(log2(expData+1));
#Read Hallmark Gene Sets & Format
hallmarkSets <- getGmt("references/h.all.v7.0.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
hallmarkSets <- geneIds(hallmarkSets);

#Run GSVA
GeneSetExprsMat <- gsva(expData, hallmarkSets, abs.ranking=F, min.sz=1, max.sz=500, parallel.sz=1, mx.diff=F);
saveRDS(GeneSetExprsMat, "results/GeneSetExpressionMatrix.RDS")
saveRDS(GeneSetExprsMat, "../../scratch/GeneSetExpressionMatrix.RDS")

