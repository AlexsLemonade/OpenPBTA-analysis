##########################################
#Purpose: Code to analyze and plot ssGSEA Results
#Author: Pichai Raman
#Date: 9/23/2019
##########################################

#Parameters for analysis
args = commandArgs(trailingOnly=TRUE) 
aovParam <- ifelse(is.na(args[1]), 0.01, as.numeric(args[1]))
padjParam <- ifelse(is.na(args[2]), 0.05, as.numeric(args[2]))
print(paste("anova P-value threshold is", aovParam))
print(paste("Tukey HSD Adjusted P-value threshold is", padjParam))

#Call libraries
library("tidyverse")
library("pheatmap")
library("viridis")

#Read in data
clinData <- read.delim("../../data/pbta-histologies.tsv", stringsAsFactors=F)
geneSetExpMat <- readRDS("../../scratch/GeneSetExpressionMatrix.RDS")

#Choose Clinical Samples corresponding to samples in Gene Set Matrix
clinData <- clinData[clinData[,"experimental_strategy"]=="RNA-Seq",]
rownames(clinData) <- clinData[,"Kids_First_Biospecimen_ID"]
clinData <- clinData[colnames(geneSetExpMat),]

#Function to split pathway names from TukeyHSD into vector of 2 elements
#i.e. Tukey HSD outputs Pathway1-Pathway2 and this will spitout c(Pathway1, Pathway2)
#Will be used so that we can tell quickly which pathways are consistently upregulated vs other pathways
getVarNames <- function(x)
{
	tmp <- strsplit(x, split="-")[[1]]
	tmp1 <- tmp[1]
	tmp2 <- tmp[2]
	return(c(tmp1, tmp2))
}

#Analyse Pathways statistically correlated with clinical variables
#Returns ANOVA P-value and Tukey HSD Results
#Requires Pathway (pathwayName) and Clinical Variable (clinVarName)
compareClinVarToPathway <- function(pathwayName=NULL, clinVarName=NULL)
{
	tmpDat <- data.frame(clinData[,clinVarName], geneSetExpMat[pathwayName,])
	colnames(tmpDat) <- c("ClinVar", "PathwayVar")
	aovOut <- aov(PathwayVar ~ ClinVar, data=tmpDat)
	posthoc <- TukeyHSD(x=aovOut, conf.level=0.95)
	posthoc <- data.frame(posthoc[[1]])
	aovOutNamed <- unlist(summary(aovOut))
	posthoc[,"F_Value"] <- aovOutNamed["F value1"]
	posthoc[,"AOV_P"] <- aovOutNamed["Pr(>F)1"]
	posthoc[,"Pathway"] <- pathwayName
	posthoc[,"ClinVar"] <- clinVarName
	posthoc <- cbind(t(sapply(rownames(posthoc), FUN=getVarNames)), posthoc)
	colnames(posthoc)[1:2] <- c("VarX", "VarY")
	return(posthoc)
}

##############################
#Analysis by Disease
##############################
diseaseType <- lapply(rownames(geneSetExpMat), FUN=compareClinVarToPathway, clinVarName="short_histology")
diseaseType <- do.call("rbind", diseaseType)

#Filter to significant entries
diseaseTypeFilt <- diseaseType %>% dplyr::filter( AOV_P < aovParam, p.adj < padjParam)

#Now make all comparisons the same direction & write out results
diseaseTypeFilt[,"DiseaseX"] <- ifelse(diseaseTypeFilt[,"diff"]>0, as.character(diseaseTypeFilt[,"VarX"]), as.character(diseaseTypeFilt[,"VarY"]))
diseaseTypeFilt[,"DiseaseY"] <- ifelse(diseaseTypeFilt[,"diff"]>0, as.character(diseaseTypeFilt[,"VarY"]), as.character(diseaseTypeFilt[,"VarX"]))
diseaseTypeFilt <- diseaseTypeFilt[-1:-2]
diseaseTypeFilt <- diseaseTypeFilt[, c("DiseaseX", "DiseaseY", colnames(diseaseTypeFilt)[1:8])]
write.table(diseaseTypeFilt, "results/DiseaseCorrelationPathway.txt", sep="\t", row.names=F)

#Now get diseases highly associated with a certain pathway
diseaseTableTmp <- data.frame(table(diseaseTypeFilt[,c("DiseaseX", "Pathway")]))
diseaseTableTmp <- diseaseTableTmp[diseaseTableTmp[,"Freq"]>(nrow(geneSetExpMat)*.25),]
diseaseTableTmp <- diseaseTableTmp[order(-diseaseTableTmp[,"Freq"]),]

#Now set up to create heatmaps
keepPathways <- unique(as.character(diseaseTableTmp[,"Pathway"]))
keepDisease <- unique(as.character(diseaseTableTmp[,"DiseaseX"]))
tmpClinData <- clinData[c("short_histology", "broad_histology")]
colnames(tmpClinData)[1] <- "Histology" 

#Function to format pathway names
formatHallmark <- function(x)
{
	x <- gsub("HALLMARK", "", x)
	x <- gsub("_", " ", x)
	x <- trimws(x)
	return(x)
}
#Format to make Pathway names more legible
rownames(geneSetExpMat) <- formatHallmark(rownames(geneSetExpMat))
keepPathways <- formatHallmark(keepPathways)

#Heatmap across all samples, only significant pathways
png("plots/HeatmapPathwaysGenes_all.png", width=1080, height=720)
pheatmap::pheatmap(geneSetExpMat[keepPathways,],
border_color="black",
color=viridis::inferno(length(geneSetExpMat) - 1), 
annotation_col=tmpClinData[c("Histology")],
show_colnames=F)
dev.off()

#Heatmap across significant histologies, and pathways
tmpClinDataTmp  <- tmpClinData[tmpClinData[,"Histology"]%in%keepDisease,]
tmpGeneSetMat <- geneSetExpMat[keepPathways,rownames(tmpClinDataTmp)]
print(paste("The Filtered Gene Set Matrix has", nrow(tmpGeneSetMat), "rows"))
print(paste("The Filtered Gene Set Matrix has", ncol(tmpGeneSetMat), "columns"))
png("plots/HeatmapPathwaysGenes_Sig.png", width=1080, height=720)
pheatmap::pheatmap(tmpGeneSetMat,
border_color="black",
color=viridis::inferno(length(geneSetExpMat) - 1), 
annotation_col=tmpClinData[c("Histology")],
show_colnames=F)
dev.off()

#Function to create boxplot for certain pathways
createBox <- function(myPathway=NULL)
{
	tmpDat <- cbind(tmpClinData, geneSetExpMat[myPathway,])
	colnames(tmpDat)[3] <- "Score"
	ggplot2::ggplot(tmpDat, aes(Histology, Score))+geom_boxplot()+theme_bw()+coord_flip()+ggtitle(paste(myPathway, "- GSVA Scores"))
	ggplot2::ggsave(paste("plots/", gsub("_", " ", myPathway), " GSVA_Boxplot.png", sep=""))
}

#Create boxplots of all pathways to show variation among histologies
sapply(keepPathways, createBox)









