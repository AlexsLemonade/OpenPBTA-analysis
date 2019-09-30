##########################################
#Purpose: Code to analyze and plot ssGSEA Results
#Author: Pichai Raman
#Date: 9/23/2019
##########################################


#Call libraries
library("tidyverse");
library("pheatmap");
library("RColorBrewer");
library("viridis")

#Read in data
clinData <- read.delim("../../data/pbta-histologies.tsv", stringsAsFactors=F);
geneSetExpMat <- readRDS("../../scratch/GeneSetExpressionMatrix.RDS")

#Choose Clinical Samples corresponding to samples in Gene Set Matrix
clinData <- clinData[clinData[,"experimental_strategy"]=="RNA-Seq",]
rownames(clinData) <- clinData[,"Kids_First_Biospecimen_ID"]
clinData <- clinData[colnames(geneSetExpMat),]

#Analyse Pathways statistically correlated with clinical variables
getVarNames <- function(x)
{
	tmp <- strsplit(x, split="-")[[1]]
	tmp1 <- tmp[1]
	tmp2 <- tmp[2];
	return(c(tmp1, tmp2))
}

compareClinVarToPathway <- function(pathwayName=NULL, clinVarName=NULL)
{
	tmpDat <- data.frame(clinData[,clinVarName], geneSetExpMat[pathwayName,]);
	colnames(tmpDat) <- c("ClinVar", "PathwayVar")
	aovOut <- aov(PathwayVar ~ ClinVar, data=tmpDat)
	posthoc <- TukeyHSD(x=aovOut, conf.level=0.95);
	posthoc <- data.frame(posthoc[[1]]);
	posthoc[,"F_Value"] <- summary(aovOut)[[1]][[4]][1];
	posthoc[,"AOV_P"] <- summary(aovOut)[[1]][[5]][1];
	posthoc[,"Pathway"] <- pathwayName;
	posthoc[,"ClinVar"] <- clinVarName;
	posthoc <- cbind(t(sapply(rownames(posthoc), FUN=getVarNames)), posthoc);
	colnames(posthoc)[1:2] <- c("VarX", "VarY")
	return(posthoc);
}

##############################
#Analysis by Disease
##############################
diseaseType <- lapply(rownames(geneSetExpMat), FUN=compareClinVarToPathway, clinVarName="short_histology");
diseaseType <- do.call("rbind", diseaseType)

#Filter to significant entries
diseaseTypeFilt <- diseaseType[diseaseType[,"AOV_P"]<0.01,]
diseaseTypeFilt <- diseaseType[diseaseType[,"p.adj"]<0.05,] 

#Now make all comparisons the same direction & write out results
diseaseTypeFilt[,"DiseaseX"] <- ifelse(diseaseTypeFilt[,"diff"]>0, as.character(diseaseTypeFilt[,"VarX"]), as.character(diseaseTypeFilt[,"VarY"]))
diseaseTypeFilt[,"DiseaseY"] <- ifelse(diseaseTypeFilt[,"diff"]>0, as.character(diseaseTypeFilt[,"VarY"]), as.character(diseaseTypeFilt[,"VarX"]))
diseaseTypeFilt <- diseaseTypeFilt[-1:-2]
diseaseTypeFilt <- diseaseTypeFilt[, c("DiseaseX", "DiseaseY", colnames(diseaseTypeFilt)[1:8])]
write.table(diseaseTypeFilt, "results/DiseaseCorrelationPathway.txt", sep="\t", row.names=F)

#Now get diseases highly associated with a certain pathway
diseaseTableTmp <- data.frame(table(diseaseTypeFilt[,c("DiseaseX", "Pathway")]));
diseaseTableTmp <- diseaseTableTmp[diseaseTableTmp[,"Freq"]>(length(rownames(geneSetExpMat))*.25),]
diseaseTableTmp <- diseaseTableTmp[order(-diseaseTableTmp[,"Freq"]),]

#Now set up to create heatmaps
keepPathways <- unique(diseaseTableTmp[,"Pathway"]);
keepDisease <- unique(diseaseTableTmp[,"DiseaseX"]);
tmpClinData <- clinData[c("short_histology", "broad_histology")];
colnames(tmpClinData)[1] <- "Histology" 
rownames(geneSetExpMat) <- gsub("HALLMARK", "", rownames(geneSetExpMat))
rownames(geneSetExpMat) <- gsub("_", " ", rownames(geneSetExpMat))
rownames(geneSetExpMat) <- trimws(rownames(geneSetExpMat))

#Let's printout the current size of the matrix
print(paste("The Gene Set Matrix has", nrow(geneSetExpMat), "rows"));
print(paste("The Gene Set Matrix has", ncol(geneSetExpMat), "columns"));

#Heatmap across all samples, only significant pathways
png("plots/HeatmapPathwaysGenes_all.png", width=1080, height=720)
pheatmap(geneSetExpMat[keepPathways,],
border_color="black",
color=inferno(length(geneSetExpMat) - 1), 
annotation_col=tmpClinData[c("Histology")],
show_colnames=F)
dev.off()

#Heatmap across significant histologies, and pathways
tmpClinDataTmp  <- tmpClinData[tmpClinData[,"Histology"]%in%keepDisease,]
tmpGeneSetMat <- geneSetExpMat[keepPathways,rownames(tmpClinDataTmp)]
png("plots/HeatmapPathwaysGenes_Sig.png", width=1080, height=720)
pheatmap(tmpGeneSetMat,
border_color="black",
color=inferno(length(geneSetExpMat) - 1), 
annotation_col=tmpClinData[c("Histology")],
show_colnames=F)
dev.off()

#Function to create boxplot for certain pathways
createBox <- function(myPathway=NULL)
{
	tmpDat <- cbind(tmpClinData, geneSetExpMat[myPathway,]);
	colnames(tmpDat)[3] <- "Score"
	ggplot(tmpDat, aes(Histology, Score))+geom_boxplot()+theme_bw()+coord_flip()+ggtitle(paste(myPathway, "- GSVA Scores"));
	ggsave(paste("plots/", gsub("_", " ", myPathway), " GSVA_Boxplot.png", sep=""))

}

#Create a few boxplots to illustrate pathways that show variation among histologies
createBox("BILE ACID METABOLISM")
createBox("E2F TARGETS")
createBox("G2M CHECKPOINT")
createBox("TGF BETA SIGNALING")
createBox("EPITHELIAL MESENCHYMAL TRANSITION")









