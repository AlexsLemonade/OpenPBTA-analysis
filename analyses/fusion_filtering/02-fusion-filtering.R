# K. S. Gaonkar 2019
# Filters standardized fusion calls to remove artifacts and false positives.
# Events such as polymerase read-throughs, mis-mapping due to gene homology, and fusions occurring in healthy normal 
# tissue require stringent filtering, making it difficult for researchers and clinicians to discern true underlying 
# oncogenic drivers of a tumor and in some cases, appropriate therapy
#
#
# Command line arguments
# 
# expressionFilter		:Integer threshold of expression for both gene in fusion  partners with FPKM<1
# junctionReadCountFilter	:Integer threshold for JunctionReadCount per fusion to remove false calls with 0 
#				supporting junction reads
# spanningFragCountFilter	:Integer threashold for (SpanningFragCount-JunctionReadCount) to remove false 
#				positives where breakpoints are towards to begining or ends of the transcript 
# reathroughFilter		:Boolean value to remove predicted read-through from caller 
# artifactFilter		:Comma separated values to remove fusions annotated as potential red flag as per
#				annotation in "annots" column (in OpenPBTA annotation is from FusionAnnotator)
# putativeDriverGeneList	:Comma separated putative driver gene list (in OpenPBTA putative gene list are 
#				TSGs,Cosmic,Oncogenic,TCGA fusion list)
# inFrame			:Boolean value to remove In-frame fusions is called by standardized fusion calls
# filterGeneList		:Comma separated gene list to capture as additional genes of interest other than 
#				putative driver list (in OpenPBTA additional filter list are Kinase, TF list)
#
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
              help="readthrough filter"),
  make_option(c("-t","--expressionFilter"), type="integer",
              help="threshold for TPM/FPKM filter"),
  make_option(c("-a","--artifactFilter"),type="character",
              help="red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator"),
  make_option(c("-j","--junctionReadCountFilter"), type="integer",
              help="threshold for JunctionReadCount"),
  make_option(c("-s","--spanningFragCountFilter"),type="integer",
              help="threshold for (SpanningFragCount - JunctionReadCount)"),
  make_option(c("-i","--readingFrameFilter"),type="character",
               help="reading frame filtering ( regex to capture inframe|frameshift|other)"),
  make_option(c("-p","--putativeDriverGeneList"),type="character",
                help="Comma separated filename for putative driver gene filter"),
  make_option(c("-f","--filterGeneList"),type="character",
                help="Comma separated filename to filter fusion with gene of interest other than putative driver list"),
  make_option(c("-R","--referenceFolder"),type="character",
                help="reference folder with required gene lists"),
  make_option(c("-o","--outputfile"),type="character",
              help="Filtered fusion calls (prefix for _putDriver.RDS and _filtFusion.RDS)")
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
putativeDriverGeneList<-opt$putativeDriverGeneList
filterGeneList<-opt$filterGeneList
readingFrameFilter<-opt$readingFrameFilter

# read standardized fusion calls
standardFusioncalls<-read.delim(standardFusionFile,stringsAsFactors=FALSE)
# to obtain geneA and geneB for gene search below
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$FusionName, pattern = '--', names = c("GeneA","GeneB")))
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$GeneA, pattern = '/', names = c("Gene1A","Gene2A")))
standardFusioncalls <- cbind(standardFusioncalls, colsplit(standardFusioncalls$GeneB, pattern = '/', names = c("Gene1B","Gene2B")))


# gather gene and fusion list from reference folder
# tab delimited gene of interest merge with filename for putativeDriver
putDriverReferenceDataTab<- tibble(filename = c(paste0(opt$referenceFolder,unlist(strsplit(putativeDriverGeneList,","))))) %>% 
  # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,read.delim,stringsAsFactor=FALSE,sep="\t")) # a new data column
putDriverReferenceDataTab<-tidyr::unnest(putDriverReferenceDataTab) %>% as.data.frame()

#gather column 1 as GeneName 2 as FusionName and 3 as source
putDriverReferenceDataTab<-data.frame("GeneName"=c(putDriverReferenceDataTab$symbol[!is.na(putDriverReferenceDataTab$symbol)],putDriverReferenceDataTab$GeneSymbol[!is.na(putDriverReferenceDataTab$GeneSymbol)],putDriverReferenceDataTab$Gene.Symbol[!is.na(putDriverReferenceDataTab$Gene.Symbol)],putDriverReferenceDataTab$Gene_A[!is.na(putDriverReferenceDataTab$Gene_A)],putDriverReferenceDataTab$Gene_B[!is.na(putDriverReferenceDataTab$Gene_B)]),"FusionName"=NA,"source"=c(putDriverReferenceDataTab$filename[!is.na(putDriverReferenceDataTab$symbol)],putDriverReferenceDataTab$filename[!is.na(putDriverReferenceDataTab$GeneSymbol)],putDriverReferenceDataTab$filename[!is.na(putDriverReferenceDataTab$Gene.Symbol)],putDriverReferenceDataTab$filename[!is.na(putDriverReferenceDataTab$Gene_A)],putDriverReferenceDataTab$filename[!is.na(putDriverReferenceDataTab$Gene_B)]),stringsAsFactors = FALSE)


# tab delimited gene of interest merge with filename for filter fusion with genes other than putative driver gene list
filtFusionrReferenceDataTab<- tibble(filename = c(paste0(opt$referenceFolder,unlist(strsplit(filterGeneList,","))))) %>% 
  # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,read.delim,stringsAsFactor=FALSE,sep="\t")) # a new data column
filtFusionrReferenceDataTab<-tidyr::unnest(filtFusionrReferenceDataTab) %>% as.data.frame()

#gather column 1 as GeneName 2 as FusionName and 3 as source
filtFusionrReferenceDataTab<-data.frame("GeneName"=c(filtFusionrReferenceDataTab$Name[!is.na(filtFusionrReferenceDataTab$Name)],filtFusionrReferenceDataTab$GeneSym[!is.na(filtFusionrReferenceDataTab$GeneSym)]),"FusionName"=NA,"source"=c(filtFusionrReferenceDataTab$filename[!is.na(filtFusionrReferenceDataTab$Name)],filtFusionrReferenceDataTab$filename[!is.na(filtFusionrReferenceDataTab$GeneSym)]),stringsAsFactors = FALSE)

fusion_filtering_Genelist<-function(standardFusioncalls=standardFusioncalls,putativeDriverGeneList=putativeDriverGeneList,filterGeneList=filterGeneList){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param putativeDriverGeneList A dataframe of genes of interest and fusions of interest ; columns 1 : GeneNames 2: FusionName 3: Source (optional)
  # @param filterGeneList A dataframe of genes of interest and fusions of interest other than putative driver list ; columns 1 : GeneNames 2: FusionName 3: Source (optional)
  # @return Standardized fusion calls filtered for genes of interest as putativeDriver and filteredFusion calls
  
  #putative driver gene and fusion list column 1 : GeneNames 2: FusionName 3: source
  putDriverFusioncalls<- standardFusioncalls %>% dplyr::filter(
    (GeneA %in% putDriverReferenceDataTab$GeneName | GeneB %in% putDriverReferenceDataTab$GeneName | FusionName %in% putDriverReferenceDataTab$FusionName)
  )
  
  #filter gene and fusion list column 1 : GeneNames 2: FusionName 3: source
  filteredFusionvcalls<- standardFusioncalls %>% dplyr::filter(
    (GeneA %in% filtFusionrReferenceDataTab$GeneName | GeneB %in% filtFusionrReferenceDataTab$GeneName | FusionName %in% filtFusionrReferenceDataTab$FusionName)
  )
  genefiltercalls<-list("filterFusion"=filteredFusionvcalls,"putDriver"=putDriverFusioncalls)
  return(genefiltercalls)
  
}


geneFiltered<-fusion_filtering_Genelist(standardFusioncalls = standardFusioncalls,putativeDriverGeneList = putDriverReferenceDataTab,filterGeneList = filtFusionrReferenceDataTab)



fusion_filtering_QC<-function(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,reathroughFilter=reathroughFilter){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param expressionMatrix A Rdata with expression data (RSEM/TPM counts) for the same cohort/GTEX
  # @param readingFramFilter A regex to capture readingframe (eg. inframe|frameshift|other)
  # @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column: "annots"
  # @param junctionReadCountFilter An integer threshold for JunctionReadCount
  # @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
  # @param expressionFilter An integer threshold for TPM/FPKM filter
  # @return Standardized fusion calls filtered to pass QC and remove calls with insufficient read-support and annotation red-flags
  
  if( reathroughFilter==TRUE & nrow(standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),]) > 0){
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

# Putative Driver
putDriverQCGenelistFiltered<-fusion_filtering_QC(standardFusioncalls = geneFiltered$putDriver,junctionReadCountFilter = junctionReadCountFilter,spanningFragCountFilter = spanningFragCountFilter,readingFrameFilter = readingFrameFilter,artifactFilter = artifactFilter,reathroughFilter = reathroughFilter)

# Filtered Fusion
filtFusionQCGenelistFiltered<-fusion_filtering_QC(standardFusioncalls = geneFiltered$filterFusion,junctionReadCountFilter = junctionReadCountFilter,spanningFragCountFilter = spanningFragCountFilter,readingFrameFilter = readingFrameFilter,artifactFilter = artifactFilter,reathroughFilter = reathroughFilter)


# load expressionMatrix RDS for expression based filtering for less than given threshold
expressionMatrix<-readRDS(expressionMatrix)
expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))


fusion_filtering_expression<-function(Sample=Sample,standardFusioncalls=standardFusioncalls,expressionFilter=expressionFilter){
  print(Sample)  
  #gene region fusions
  standardFusioncallsGeneRegion<-standardFusioncalls[is.na(standardFusioncalls$Gene2A) & is.na(standardFusioncalls$Gene2B) ,]  
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
  standardFusioncalls <- unique(standardFusioncalls[,c('LeftBreakpoint','RightBreakpoint','FusionName' , 'Sample' , 'Caller' ,'Fusion_Type' , 'JunctionReadCount' , 'SpanningFragCount' , 'Confidence' ,'annots')])
  
  return(standardFusioncalls)
  
}

#get sample list
sampleList<-unique(c(putDriverQCGenelistFiltered$Sample,filtFusionQCGenelistFiltered$Sample))

# run per sample and merge putative Driver fusions
standardFusioncallsPerSample<-lapply(as.character(sampleList),function(x) fusion_filtering_expression(Sample = x,expressionFilter = 1,standardFusioncalls = putDriverQCGenelistFiltered))
putDriverstandardFusioncalls <- Reduce(function(x, y) rbind (x,y), standardFusioncallsPerSample)
saveRDS(putDriverstandardFusioncalls,paste0(opt$outputfile,"_putDriver.RDS"))

# run per sample and merge putative Driver fusions
standardFusioncallsPerSample<-lapply(as.character(sampleList),function(x) fusion_filtering_expression(Sample = x,expressionFilter = 1,standardFusioncalls = filtFusionQCGenelistFiltered))
filtFusionstandardFusioncalls <- Reduce(function(x, y) rbind (x,y), standardFusioncallsPerSample)
saveRDS(filtFusionstandardFusioncalls,paste0(opt$outputfile,"_filtFusion.RDS"))


