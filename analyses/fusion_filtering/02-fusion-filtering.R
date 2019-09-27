# K. S. Gaonkar 2019
# Filters standardized fusion calls to remove artifacts and false positives.
# Events such as polymerase read-throughs, mis-mapping due to gene homology, and fusions occurring in healthy normal 
# tissue require stringent filtering, making it difficult for researchers and clinicians to discern true underlying 
# oncogenic drivers of a tumor and in some cases, appropriate therapy
#
#
# Command line arguments
#
# standardFusionFile      :Standardized fusion calls from [STAR|ARRIBA] from 01-fusion-standardization.R 
#				                   or user input files in the correct format 
# expressionFilter        :Integer threshold of expression for both gene in fusion  partners with FPKM<1
# junctionReadCountFilter	:Integer threshold for JunctionReadCount per fusion to remove false calls with 0 
#				                   supporting junction reads
# spanningFragCountFilter	:Integer threashold for (SpanningFragCount-JunctionReadCount) to remove false 
#				                   positives where breakpoints are towards to begining or ends of the transcript 
# reathroughFilter		    :Boolean value to remove predicted read-through from caller 
# artifactFilter		      :Comma separated values to remove fusions annotated as potential red flag as per
#				                   annotation in "annots" column (in OpenPBTA annotation is from FusionAnnotator)
# readingFrameFilter		  :Reading frame to keep in final set of QC fusion calls ( regex to capture inframe|frameshift|other)
# referenceFolder		      :Path to folder with all reference gene list and fusion file list with the following files 
#                          genelistreference.txt A dataframe of genes of interest ; columns 1 : GeneNames 2: Source file 3: Type
#                          fusionreference.txt A dataframe of fusion of interest ; columns 1 : FusionName 2: Source file 3: Type
# outputfile			        :Filename prefix for QC filtered and gene of interest annotated fusion calls (prefix for _QC_expression_filtered_annotated.RDS)")
#

  
				
				
				

library("optparse")
library("dplyr")
library("tidyr")
library("purrr")
library("reshape2")

option_list <- list(
  make_option(c("-S", "--standardFusionFile"),type="character",
              help="Standardized fusion calls from (.TSV) "),
  make_option(c("-e","--expressionMatrix"),type="character",
              help="Matrix of expression for samples in standardFusion file (.RDS) "),
  make_option(c("-r", "--reathroughFilter"), type="character",
              help="Boolean to filter readthrough",default=TRUE),
  make_option(c("-t","--expressionFilter"), type="integer",
              help="threshold for TPM/FPKM filter",default=1),
  make_option(c("-a","--artifactFilter"),type="character",
              help="red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator"),
  make_option(c("-j","--junctionReadCountFilter"), type="integer",
              help="threshold for JunctionReadCount"),
  make_option(c("-s","--spanningFragCountFilter"),type="integer",
              help="threshold for (SpanningFragCount - JunctionReadCount)"),
  make_option(c("-i","--readingFrameFilter"),type="character",
               help="reading frame filtering ( regex to capture inframe|frameshift|other)"),
  make_option(c("-R","--referenceFolder"),type="character",
                help="reference folder with required gene lists"),
  make_option(c("-o","--outputfile"),type="character",
              help="Filtered fusion calls (prefix for _QC_expression_filtered_annotated.RDS)")
)


# get command line options, if help option encountered print help and exit,

opt <- parse_args(OptionParser(option_list=option_list))

standardFusionFile<-opt$standardFusionFile
expressionMatrix<-opt$expressionMatrix
reathroughFilter<-opt$reathroughFilter
expressionFilter<-opt$expressionFilter
artifactFilter<-opt$artifactFilter
junctionReadCountFilter<-opt$junctionReadCountFilter
spanningFragCountFilter<-opt$spanningFragCountFilter
referenceFolder<-opt$referenceFolder
readingFrameFilter<-opt$readingFrameFilter

# read standardized fusion calls
standardFusioncalls<-read.delim(standardFusionFile,stringsAsFactors=FALSE)
# to obtain geneA and geneB for gene search below
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$FusionName, pattern = '--', names = c("GeneA","GeneB")))
# Intergenic fusion will have Gene1A,Gene2A,Gene1B,Gene2B
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$GeneA, pattern = '/', names = c("Gene1A","Gene2A")))
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$GeneB, pattern = '/', names = c("Gene1B","Gene2B")))

#formatting dataframe for filtering
standardFusioncalls<-standardFusioncalls %>% 
  #For annotation later
  mutate (Gene1A_anno=NA,Gene2A_anno=NA,Gene1B_anno=NA,Gene2B_anno=NA) %>% 
  #remove distance to fusion breakpoint from gene names in intergenic fusions
  mutate(Gene1A=gsub("[(].*","",standardFusioncalls$Gene1A),
         Gene2A=gsub("[(].*","",standardFusioncalls$Gene2A),
         Gene1B=gsub("[(].*","",standardFusioncalls$Gene1B),
         Gene2B=gsub("[(].*","",standardFusioncalls$Gene2B))



#############artifactFiltering#############


fusion_filtering_QC<-function(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,reathroughFilter=reathroughFilter){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param expressionMatrix A Rdata with expression data (RSEM/TPM counts) for the same cohort/GTEX
  # @param readingFramFilter A regex to capture readingframe (eg. inframe|frameshift|other)
  # @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column: "annots"
  # @param junctionReadCountFilter An integer threshold for JunctionReadCount
  # @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
  # @param expressionFilter An integer threshold for TPM/FPKM filter
  # @return Standardized fusion calls filtered to pass QC and remove calls with insufficient read-support and annotation red-flags
  
  if( reathroughFilter!=FALSE & nrow(standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),]) > 0){
    # Gather read throughs from standardized fusion calls  
    rts <- standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),"FusionName"]
    # Reverse of read throughs to capture 
    rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
    # Combine read through and reverse fusion genes 
    rts <- unique(c(rts, rts.rev))
    standardFusioncalls <- standardFusioncalls[-which(standardFusioncalls$FusionName %in% rts),]
  }

  if( !missing(readingFrameFilter ) ){
    standardFusioncalls <- standardFusioncalls[grep(readingFrameFilter,standardFusioncalls$Fusion_Type),]
  }
  
  if( !missing(artifactFilter) & nrow(standardFusioncalls[grep(artifactFilter,standardFusioncalls$annots),])>0 ) {
    standardFusioncalls <- standardFusioncalls[-grep(artifactFilter,standardFusioncalls$annots),]
  }

  if( !missing(junctionReadCountFilter ) ){
    standardFusioncalls <- standardFusioncalls[which(standardFusioncalls$JunctionReadCount>junctionReadCountFilter),]
  }

  if( !missing(spanningFragCountFilter ) ){
    # false positive calls either at the start or end of a transcript will be supported by an uneven majority of support from spanning fragments compared to junction reads
    # to remove these calls we are implementing this condition below
    standardFusioncalls <- standardFusioncalls[which( (standardFusioncalls$SpanningFragCount- standardFusioncalls$JunctionReadCount ) < spanningFragCountFilter),]
  }
  
  return(standardFusioncalls)
}

# QC filter: artifact and read support
QCFiltered<-fusion_filtering_QC(standardFusioncalls = standardFusioncalls,junctionReadCountFilter = junctionReadCountFilter,spanningFragCountFilter = spanningFragCountFilter,readingFrameFilter = readingFrameFilter,artifactFilter = artifactFilter,reathroughFilter = reathroughFilter)


# load expressionMatrix RDS for expression based filtering for less than given threshold
expressionMatrix<-readRDS(expressionMatrix)
expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))


fusion_filtering_expression<-function(Sample=Sample,standardFusioncalls=standardFusioncalls,expressionFilter=expressionFilter){
  print(Sample)  
  #gene region fusions
  standardFusioncallsGeneRegion<-standardFusioncalls[which(standardFusioncalls$Gene2A=="" & standardFusioncalls$Gene2B=="") ,]  
  #add column with expression for geneA
  standardFusioncallsGeneA<-merge(standardFusioncallsGeneRegion[which(standardFusioncallsGeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="GeneA",by.y="GeneSymbol")
  #add column with expression for geneB
  standardFusioncallsGeneB<-merge(standardFusioncallsGeneRegion[which(standardFusioncallsGeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="GeneB",by.y="GeneSymbol")
  standardFusioncallsGeneRegion<-rbind(standardFusioncallsGeneA,standardFusioncallsGeneB)
  
  
  
  #intergenic gene fusions
  standardFusioncallsIntergeneRegion<-standardFusioncalls[grep("/",standardFusioncalls$FusionName),]  
  if(nrow(standardFusioncallsIntergeneRegion)>0){
    #add column with expression for gene1A
    standardFusioncallsGene1A<-merge(standardFusioncallsIntergeneRegion[which(standardFusioncallsIntergeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="Gene1A",by.y="GeneSymbol")
    #add column with expression for gene1B
    standardFusioncallsGene1B<-merge(standardFusioncallsIntergeneRegion[which(standardFusioncallsIntergeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="Gene1B",by.y="GeneSymbol")
    #add column with expression for gene2A
    standardFusioncallsGene2A<-merge(standardFusioncallsIntergeneRegion[which(standardFusioncallsIntergeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="Gene2A",by.y="GeneSymbol")
    #add column with expression for gene2B
    standardFusioncallsGene2B<-merge(standardFusioncallsIntergeneRegion[which(standardFusioncallsIntergeneRegion$Sample==Sample),],expressionMatrix[,c("GeneSymbol",Sample)],by.x="Gene2B",by.y="GeneSymbol")
    standardFusioncallsIntergeneRegion<-rbind(standardFusioncallsGene1A,standardFusioncallsGene2A,standardFusioncallsGene1B,standardFusioncallsGene2B)
  }
  
  #merge gene and intergenic fusions
  standardFusioncalls<-rbind(standardFusioncallsGeneRegion,standardFusioncallsIntergeneRegion)
  
  #filter < expressionthreshold
  if( !missing(expressionFilter ) ){
    standardFusioncalls <- standardFusioncalls[which(as.integer(standardFusioncalls[,Sample]) > expressionFilter),]
  }
  
  #get unique fusion calls 
  standardFusioncalls <- unique(standardFusioncalls[,c('LeftBreakpoint','RightBreakpoint','FusionName' , 'Sample' , 'Caller' ,'Fusion_Type' , 'JunctionReadCount' , 'SpanningFragCount' , 'Confidence' ,'annots','Gene1A','Gene2A','Gene1B','Gene2B')])
  
  return(standardFusioncalls)
  
}

#get sample list
sampleList<-unique(QCFiltered$Sample)

# run per sample on standard fusion calls
standardFusioncallsPerSample<-lapply(as.character(sampleList),function(x) fusion_filtering_expression(Sample = x,expressionFilter = expressionFilter,standardFusioncalls = QCFiltered))

QCExpFilteredstandardFusioncalls <- Reduce(function(x, y) rbind (x,y), standardFusioncallsPerSample)
saveRDS(QCExpFilteredstandardFusioncalls,paste0(opt$outputfile,"_QC_expression_filtered.RDS"))

#####################################


###############annotation############
# column 1 as GeneName 2 source file 3 Type; collapse to summarize type
geneListReferenceDataTab<-read.delim(paste0(referenceFolder,"genelistreference.txt"),stringsAsFactors = FALSE)
geneListReferenceDataTab<-geneListReferenceDataTab %>% group_by(Gene_Symbol) %>% mutate(type = toString(type)) %>% as.data.frame()
geneListReferenceDataTab<-unique(geneListReferenceDataTab[,c("Gene_Symbol","type")])

# column 1 as FusionName 2 source file 3 Type; collapse to summarize type
fusionReferenceDataTab<-read.delim(paste0(referenceFolder,"fusionreference.txt"),stringsAsFactors = FALSE)
fusionReferenceDataTab<-unique(fusionReferenceDataTab[,c("FusionName","type")])

fusion_annotation_list<-function(standardFusioncalls=standardFusioncalls,geneListReference=geneListReference,fusionReference=fusionReference){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized and QC filtered 
  # @param geneListReference A dataframe of genes of interest ; columns 1 : GeneNames 2: Source file 3: Type
  # @param fusionReference A dataframe of fusion of interest ; columns 1 : FusionName 2: Source file 3: Type
  # @return Standardized and QC filtered fusion calls annotated with gene list and fusion calls from reference lists
  
#Gene annotation
standardFusioncallsGOI<-standardFusioncalls %>% 
           mutate(Gene1A_anno=lapply(standardFusioncalls$Gene1A,function(x) geneListReferenceDataTab[which(geneListReferenceDataTab$Gene_Symbol==x),"type"]),
           Gene2A_anno=lapply(standardFusioncalls$Gene2A,function(x) geneListReferenceDataTab[which(geneListReferenceDataTab$Gene_Symbol==x),"type"]),
           Gene1B_anno=lapply(standardFusioncalls$Gene1B,function(x) geneListReferenceDataTab[which(geneListReferenceDataTab$Gene_Symbol==x),"type"]),
           Gene2B_anno=lapply(standardFusioncalls$Gene2B,function(x) geneListReferenceDataTab[which(geneListReferenceDataTab$Gene_Symbol==x),"type"]))
  
#Fusion annotation  
standardFusioncallsFOI<-standardFusioncalls %>% 
  mutate(Fusion_anno=lapply(standardFusioncalls$FusionName,function(x) fusionReferenceDataTab[which(fusionReferenceDataTab$FusionName==x),"type"]))

standardFusioncalls<-list("geneAnnotion"=standardFusioncallsGOI,"fusionAnnottaion"=standardFusioncallsFOI)

return(standardFusioncalls)
}


standardFusioncallsAnnotated<-fusion_annotation_list(standardFusioncalls = standardFusioncalls,geneListReference=geneListReference,fusionReference=fusionReference)

saveRDS(standardFusioncallsAnnotated,paste0(opt$outputfile,"_QC_expression_filtered_annotated.RDS"))


######################################

