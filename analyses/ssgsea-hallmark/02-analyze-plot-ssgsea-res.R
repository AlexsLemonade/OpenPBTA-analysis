##########################################
#Purpose: Code to analyze and plot ssGSEA Results
#Author: Pichai Raman
#Date: 9/23/2019
##########################################

#Call libraries
library("tidyverse")
library("pheatmap")
library("viridis")
library("optparse")

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-a", "--anova_pvalue"),
    type = "double",
    default = 0.01,
    help = "ANOVA P-value threshold",
  ),
  optparse::make_option(
    c("-t", "--tukey_pvalue"),
    type = "double",
    default = 0.05,
    help = "Tukey HSD Adjusted P-value threshold"
  ),
  optparse::make_option(
    c("-p", "--perc_keep"),
    type = "double",
    default = .25,
    help = "Percentage of diseases a particular disease has a significantly higher score than for a specific pathway"
  )
)

#Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#Pass arguments to variables
aovParam <- opt$anova_pvalue
padjParam <- opt$tukey_pvalue
percKeep <- opt$perc_keep


#Read in clinical data and pathway level data generated from 01-run-ssgsea.R
clinData <- readr::read_tsv("../../data/pbta-histologies.tsv", guess_max = 10000)
geneSetExpMat <- readRDS("../../scratch/GeneSetExpressionMatrix.RDS")

#Choose Clinical Samples corresponding to samples in Gene Set Matrix
#This is to makes sure that clinical data and the pathway DF have the same samples
#and that they are ordered the same way
clinData <- clinData[clinData[,"experimental_strategy"]=="RNA-Seq",]
rownames(clinData) <- clinData[,"Kids_First_Biospecimen_ID"]
clinData <- clinData[colnames(geneSetExpMat),]

#Function to split variable names from TukeyHSD into vector of 2 elements
#i.e. Tukey HSD outputs variable1-variable2 and this will spitout c(variable1, variable2)
#This will make for easier filtering and/or summarization on one variable
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
	colnames(tmpDat) <- c("ClinVar", "PathwayVar") 	#Need to set this since pathwayName and clinVarName will change which is not good for aov function
	aovOut <- aov(PathwayVar ~ ClinVar, data=tmpDat)
	posthoc <- TukeyHSD(x=aovOut, conf.level=0.95)
	posthoc <- data.frame(posthoc[[1]])
	aovOutNamed <- unlist(summary(aovOut))
	posthoc[,"F_Value"] <- aovOutNamed["F value1"]
	posthoc[,"AOV_P"] <- aovOutNamed["Pr(>F)1"]
	posthoc[,"Pathway"] <- pathwayName
	posthoc[,"ClinVar"] <- clinVarName
	posthoc <- cbind(t(sapply(rownames(posthoc), FUN=getVarNames)), posthoc) #Pulls out clinical variable names 
	colnames(posthoc)[1:2] <- c("VarX", "VarY")
	return(posthoc)
}


#Analysis by short_histology. This is only clinical variable we are examining but the function compareClinVarToPathway
#is generic enough that any variable in clinical data can be used. This will be useful if we decide to use for instance molecular
#subtypes instead of histology for this analysis. So the lapply below will for each pathway
#run an ANOVA across all histologies i.e. looking to see if any particular histology has a significantly greater score
#than any other histology and then follow that with a post-hoc test. The comparison is not done between pathways. 
diseaseType <- lapply(rownames(geneSetExpMat), FUN=compareClinVarToPathway, clinVarName="short_histology")
diseaseType <- do.call("rbind", diseaseType) 

#Filter to significant entries by F and Tukey and write out 
diseaseTypeFilt <- diseaseType %>% dplyr::filter( AOV_P <= aovParam, p.adj <= padjParam)
colnames(diseaseTypeFilt)[1:2] <- c("DiseaseX", "DiseaseY") #Rename columns to match clinical variable type
write.table(diseaseTypeFilt, "results/DiseaseCorrelationPathway.txt", sep="\t", row.names=F)

#We would like to filter to cases where a histology is significantly higher than other histologies (and the corresponding pathways)
#This is for the heatmap so we may focus a bit on the more interesting results
#The issue here is that Tukey will write out the histologies in some particular order
#If we are interested in finding if a Histology A is significantly upregulated compared to at least X% of
#other histologies this may present a problem because Tukey may write out the following
#HistologyA-HistologyB ~ Fold-change = 1
#HistologyA-HistologyC ~ Fold-change = 1
#HistologyD-HistologyA ~ Fold-change = -1
#if the first histology is in column 1 and the second is in column 2 and the FC is in column 3 of a data frame we would be required
#to check the FC in order to determine if Histology Y is in fact > or < another histology. Hence a bit of logic below to get counts. 
diseaseTypeFilt[,"DiseaseWithHigherAvgScore"] <- ifelse(diseaseTypeFilt[,"diff"]>0, as.character(diseaseTypeFilt[,"DiseaseX"]), as.character(diseaseTypeFilt[,"DiseaseY"]))
diseaseTableTmp <- data.frame(table(diseaseTypeFilt[,c("DiseaseWithHigherAvgScore", "Pathway")])) 

#Now get diseases highly associated with a certain pathway based on a cutoff, which is based on the number of histologies in general
diseaseTableTmp <- diseaseTableTmp[diseaseTableTmp[,"Freq"]>(n_distinct(clinData[,"short_histology"])*percKeep),]
diseaseTableTmp <- diseaseTableTmp[order(-diseaseTableTmp[,"Freq"]),]

#Now set up to create heatmaps
keepPathways <- unique(as.character(diseaseTableTmp[,"Pathway"]))
keepDisease <- unique(as.character(diseaseTableTmp[,"DiseaseWithHigherAvgScore"]))
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

#Pathways to keep
print(paste("Number of pathways is", dim(geneSetExpMat[keepPathways,])))
print(keepPathways)

#Heatmap across all samples, only significant pathways
png("plots/HeatmapPathwaysGenes_all.png", width=1080, height=720)
	pheatmap::pheatmap(geneSetExpMat[keepPathways,],
	border_color="black",
	color=viridis::inferno(length(geneSetExpMat) - 1), 
	annotation_col=tmpClinData[c("Histology")],
	show_colnames=F)
dev.off()

#Heatmap across significant histologies and significant pathways
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
	ggplot2::ggplot(tmpDat, aes(Histology, geneSetExpMat[myPathway, ]))+
	ggplot2::geom_boxplot()+
	ggplot2::theme_bw()+
	ggplot2::coord_flip()+
	ggplot2::ggtitle(paste(myPathway, "- GSVA Scores"))+
	ggplot2::xlab("Score")
	ggplot2::ggsave(paste("plots/", gsub("_", " ", myPathway), " GSVA_Boxplot.png", sep=""))
}

#Create boxplots of all pathways to show variation among histologies
sapply(keepPathways, createBox)









