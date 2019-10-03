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
  make_option(c("-r", "--reathroughFilter"), type="character",action = "store_true",
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
  #remove distance to fusion breakpoint from gene names in intergenic fusions
  mutate(Gene1A=gsub("[(].*","",standardFusioncalls$Gene1A),
         Gene2A=gsub("[(].*","",standardFusioncalls$Gene2A),
         Gene1B=gsub("[(].*","",standardFusioncalls$Gene1B),
         Gene2B=gsub("[(].*","",standardFusioncalls$Gene2B))



#############artifactFiltering#############


fusion_filtering_QC<-function(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,reathroughFilter=reathroughFilter){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param expressionMatrix A rds with expression data (RSEM/TPM counts) for the same cohort/GTEX
  # @param readingFramFilter A regex to capture readingframe (eg. inframe|frameshift|other)
  # @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column: "annots"
  # @param junctionReadCountFilter An integer threshold for JunctionReadCount
  # @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
  # @param expressionFilter An integer threshold for TPM/FPKM filter
  # @return Standardized fusion calls filtered to pass QC and remove calls with insufficient read-support and annotation red-flags
  
  if( reathroughFilter & any(grepl('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE))){
    # Gather read throughs from standardized fusion calls  
    rts <- standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),"FusionName"]
    # Reverse of read throughs to capture 
    rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
    # Combine read through and reverse fusion genes 
    rts <- unique(c(rts, rts.rev))
    standardFusioncalls <- standardFusioncalls[-which(standardFusioncalls$FusionName %in% rts),]
  }

  if( !missing(readingFrameFilter ) ){
    #TODO check for values
    standardFusioncalls <- standardFusioncalls[grep(readingFrameFilter,standardFusioncalls$Fusion_Type),]
  }
  
  if( !missing(artifactFilter) & any(grepl(artifactFilter,standardFusioncalls$annots)) ) {
    #TODO check for values
    standardFusioncalls <- standardFusioncalls[-grep(artifactFilter,standardFusioncalls$annots),]
  }

  if( !missing(junctionReadCountFilter ) ){
    standardFusioncalls <- standardFusioncalls[which(standardFusioncalls$JunctionReadCount>=junctionReadCountFilter),]
  }

  if( !missing(spanningFragCountFilter ) ){
    # false positive calls either at the start or end of a transcript will be supported by an uneven majority of support from spanning fragments compared to junction reads
    # to remove these calls we are implementing this condition below
    standardFusioncalls <- standardFusioncalls[which( (standardFusioncalls$SpanningFragCount - standardFusioncalls$JunctionReadCount ) <= spanningFragCountFilter),]
  }
  
  return(standardFusioncalls)
}

# QC filter: artifact and read support
QCFiltered<-fusion_filtering_QC(standardFusioncalls = standardFusioncalls,junctionReadCountFilter = junctionReadCountFilter,spanningFragCountFilter = spanningFragCountFilter,readingFrameFilter = readingFrameFilter,artifactFilter = artifactFilter,reathroughFilter = reathroughFilter)

############################################

########### Expression filtering ###########

# load expressionMatrix RDS for expression based filtering for less than given threshold
expressionMatrix<-readRDS(expressionMatrix)
expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))


# The idea from @jaclyn-taroni 
# Generate two data frames that keep track of all gene symbols involved for each sample-fusion name pair and contain a sample's expression value for a gene symbol in long format. Use these to  filter all the available fusions to just the ones with either gene from the fusion pair to have expression above the threshold.

fusion_sample_gene_df <- QCFiltered %>%
  # We want to keep track of the gene symbols for each sample-fusion pair
  dplyr::select(Sample, FusionName, Gene1A, Gene1B, Gene2A, Gene2B) %>%
  # We want a single column that contains the gene symbols
  tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                key = gene_position, value = GeneSymbol) %>%
  # Get rid of the Gene1A, Gene1B, Gene2A, Gene2B information              
  dplyr::select(-gene_position) %>%
  # Remove columns without gene symbols
  dplyr::filter(GeneSymbol != "") %>%
  # This is for illustrations sake, only
  dplyr::arrange(Sample, FusionName) %>%
  # Retain only distinct rows
  dplyr::distinct()


expression_long_df <- expressionMatrix %>%
  # Keep the gene symbols and the samples themselves
  dplyr::select(GeneSymbol, dplyr::starts_with("BS_")) %>%
  # Get the data into long format
  reshape2::melt(variable.name = "Sample",
                 value.name = "expression_value") %>%
  # Remove rows with expression that is too low
  dplyr::filter(expression_value > expressionFilter)


expression_filtered_fusions <- fusion_sample_gene_df %>%
  # join the filtered expression values to the data frame keeping track of symbols
  # for each sample-fusion name pair
  dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol"))  %>%
  # for each sample-fusion name pair, are all genes under the expression threshold?
  # keep track in `all_low_expression` column
  dplyr::group_by(FusionName, Sample) %>%
  dplyr::mutate(all_low_expression = all(is.na(expression_value))) %>%
  # only keep the rows that *don't* have all below threshold
  dplyr::filter(!all_low_expression) %>%
  # we only need the FusionName and Sample to filter the entire fusion data.frame
  dplyr::select(FusionName, Sample) %>%
  # unique FusionName-Sample rows
  dplyr::distinct() %>%
  # use this to filter the QC filtered fusion data frame 
  dplyr::inner_join(QCFiltered, by = c("FusionName", "Sample")) %>% 
  # retain the same columns as merge method
  dplyr::select(c('LeftBreakpoint','RightBreakpoint','FusionName' , 'Sample' , 
                  'Caller' ,'Fusion_Type' , 'JunctionReadCount' ,'SpanningFragCount' , 
                  'Confidence' ,'annots','Gene1A','Gene2A','Gene1B','Gene2B'))



###############annotation###############
# column 1 as GeneName 2 source file 3 Type; collapse to summarize type
geneListReferenceDataTab<-read.delim(paste0(referenceFolder,"genelistreference.txt"),stringsAsFactors = FALSE)
geneListReferenceDataTab<-geneListReferenceDataTab %>% dplyr::group_by(Gene_Symbol) %>% dplyr::mutate(type = toString(type)) %>% as.data.frame()
geneListReferenceDataTab<-unique(geneListReferenceDataTab[,c("Gene_Symbol","type")])

# column 1 as FusionName 2 source file 3 Type; collapse to summarize type
fusionReferenceDataTab<-read.delim(paste0(referenceFolder,"fusionreference.txt"),stringsAsFactors = FALSE)
fusionReferenceDataTab<-unique(fusionReferenceDataTab[,c("FusionName","type")])


# Annotate QC filtered fusion calls
annotated_filtered_fusions<-expression_filtered_fusions %>% 
  # annotate Gene1A
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene1A"="Gene_Symbol")) %>% dplyr::rename(Gene1A_anno=type) %>% 
  # annotate Gene1B
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene1B"="Gene_Symbol")) %>% dplyr::rename(Gene1B_anno=type) %>% 
  # annotate Gene2A
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene2A"="Gene_Symbol")) %>% dplyr::rename(Gene2A_anno=type) %>% 
  # annotate Gene2B
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene2B"="Gene_Symbol")) %>% dplyr::rename(Gene2B_anno=type) %>% 
  # annotate FusionName
  dplyr::left_join(fusionReferenceDataTab,by=c("FusionName"="FusionName")) %>% dplyr::rename(Fusion_anno=type) %>%
  as.data.frame()

saveRDS(annotated_filtered_fusions,paste0(opt$outputfile,"_QC_expression_filtered_annotated.RDS"))

# Annotate standardized fusion calls for project specific gene list filtering
annotated_unfiltered_fusions<-standardFusioncalls %>% 
  # annotate Gene1A
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene1A"="Gene_Symbol")) %>% dplyr::rename(Gene1A_anno=type) %>% 
  # annotate Gene1B
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene1B"="Gene_Symbol")) %>% dplyr::rename(Gene1B_anno=type) %>% 
  # annotate Gene2A
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene2A"="Gene_Symbol")) %>% dplyr::rename(Gene2A_anno=type) %>% 
  # annotate Gene2B
  dplyr::left_join(geneListReferenceDataTab,by=c("Gene2B"="Gene_Symbol")) %>% dplyr::rename(Gene2B_anno=type) %>% 
  # annotate FusionName
  dplyr::left_join(fusionReferenceDataTab,by=c("FusionName"="FusionName")) %>% dplyr::rename(Fusion_anno=type) %>%
  as.data.frame()


saveRDS(annotated_unfiltered_fusions,paste0(opt$outputfile,"_unfiltered_annotated.RDS"))


############################################################################################

############################################################################################



