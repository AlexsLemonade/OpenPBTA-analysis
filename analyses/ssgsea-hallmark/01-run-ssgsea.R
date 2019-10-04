##########################################
#Purpose: Code to run ssGSEA analysis
#Author: Pichai Raman
#Date: 9/23/2019
##########################################


#Call libraries
library("GSVA")
library("GSEABase")
library("stringr")

#Read in data 
expData <- readRDS("../../data/pbta-gene-counts-rsem-expected_count.stranded.rds")

#Format RNA-Seq data and get one gene symbol per row as there is a many-to-one mapping
#of gene symbols and Ensembl ID's
expData[,"max"] <- apply(expData[-1], FUN=max, MARGIN=1) #Get max value per row
expData[,"Gene"] <- stringr::word(as.character(expData[,1]), 2, sep = "_") #Strip out gene symbol
expData <- expData[order(-expData[,"max"]),] #Order by Max Value
expData <- expData[!duplicated(expData[,"Gene"]),] #Remove duplicates - i.e we are choosing the Ensembl ID with the max value as our representative gene symbol
rownames(expData) <- expData[,"Gene"] #Set Rownames to gene symbols instead of genesymbol_Ensembl
expData <- expData[expData[,"max"]>0,] #Also get rid of any rows with an FPKM of 0 for all samples
expData <- dplyr::select(expData, dplyr::starts_with("BS_")) #Now we have on symbol per row in the matrix

#Log-Tranform for input into GSVA
expData <- as.matrix(log2(expData+1))

#Read Hallmark Gene Sets (function getGMT from GSEABase) and use geneIDs function for input into GSVA
hallmarkSets <- GSEABase::getGmt("references/h.all.v7.0.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
hallmarkSets <- GSEABase::geneIds(hallmarkSets)

#Run GSVA to get Hallmark pathway level scores for each sample
GeneSetExprsMat <- GSVA::gsva(expData, hallmarkSets, abs.ranking=F, min.sz=1, max.sz=500, parallel.sz=1, mx.diff=F)
saveRDS(GeneSetExprsMat, "results/GeneSetExpressionMatrix.RDS") #Write to results 
saveRDS(GeneSetExprsMat, "../../scratch/GeneSetExpressionMatrix.RDS") #Write to scratch for use in next function

