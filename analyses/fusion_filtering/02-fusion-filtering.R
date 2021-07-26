# K. S. Gaonkar 2019
# Filters standardized fusion calls to remove artifacts and false positives.
# Events such as polymerase read-throughs, mis-mapping due to gene homology, and fusions occurring in healthy normal
# tissue require stringent filtering, making it difficult for researchers and clinicians to discern true underlying
# oncogenic drivers of a tumor and in some cases, appropriate therapy

# Example run:
# Rscript analyses/fusion_filtering/02-fusion-filtering.R -S scratch/arriba.tsv --expressionMatrix data/gene-expression-rsem-tpm-collapsed.rds --readthroughFilter --artifactFilter "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap" --junctionReadCountFilter 1 --spanningFragCountFilter 100 --readingFrameFilter "in-frame|frameshift|other" --referenceFolder analyses/fusion_filtering/references/ --outputFile scratch/standardFusionExp -t 1
#
# Command line arguments
#
# standardFusionFiles  :Standardized fusion calls from STARFusion and Arriba from 01-fusion-standardization.R as comma separated file names
# expressionFilter        :Integer threshold of expression for both gene in fusion  partners with TPM<1
# junctionReadCountFilter	:Integer threshold for JunctionReadCount per fusion to remove false calls with 0
#				                   supporting junction reads
# spanningFragCountFilter	:Integer threashold for (SpanningFragCount-JunctionReadCount) to remove false
#				                   positives where breakpoints are towards to begining or ends of the transcript
# readthroughFilter		    :Boolean value to remove predicted read-through from caller
# artifactFilter		      :Comma separated values to remove fusions annotated as potential red flag as per
#				                   annotation in "annots" column (in OpenPBTA annotation is from FusionAnnotator)
# readingFrameFilter		  :Reading frame to keep in final set of QC fusion calls ( regex to capture inframe|frameshift|other)
# referenceFolder		      :Path to folder with all reference gene list and fusion file list with the following files
#                          genelistreference.txt A dataframe of genes of interest ; columns 1 : GeneNames 2: Source file 3: Type
#                          fusionreference.txt A dataframe of fusion of interest ; columns 1 : FusionName 2: Source file 3: Type
# outputFile			        :Filename prefix for QC filtered and gene of interest annotated fusion calls (prefix for _QC_expression_filtered_annotated.RDS)")
#



library("optparse")
library("reshape2")
library("tidyverse")
library("qdapRegex")

option_list <- list(
  make_option(c("-S", "--standardFusionFiles"),type="character",
              help="Standardized fusion calls from STARFusion (.TSV), Arriba (.TSV) "),
  make_option(c("-y", "--cohortInterest"),type="character",
              help="cohort of interest for the filtering"),
  make_option(c("-e","--expressionMatrix"),type="character",
              help="Matrix of expression for samples in standardFusion file (.RDS) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-r", "--readthroughFilter"),action = "store_true",
              help="Boolean to filter readthrough",default=FALSE),
  make_option(c("-t","--expressionFilter"), type="integer",
              help="threshold for TPM/FPKM filter",default=1),
  make_option(c("-a","--artifactFilter"),type="character",
              help="red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator"),
  make_option(c("-j","--junctionReadCountFilter"), type="integer",
              help="threshold for junctionReadCount",default=1),
  make_option(c("-s","--spanningFragCountFilter"),type="integer",
              help="threshold for  (SpanningFragCount - JunctionReadCount)",default=100),
  make_option(c("-i","--readingFrameFilter"),type="character",
               help="reading frame filtering ( regex to capture inframe|frameshift|other)"),
  make_option(c("-R","--referenceFolder"),type="character",
                help="reference folder with required gene lists"),
  make_option(c("-o","--outputFile"),type="character",
              help="Filtered fusion calls (prefix for _QC_expression_filtered_annotated.RDS)")
)


# get command line options, if help option encountered print help and exit,

opt <- parse_args(OptionParser(option_list=option_list))

# TO-DO
# multiple opt values for each caller and output name from input fusion calls and expression Matrix?
standardFusionFiles<-unlist(strsplit(opt$standardFusionFiles,","))
cohortInterest<-unlist(strsplit(opt$cohortInterest,","))
expressionMatrix<-opt$expressionMatrix
readthroughFilter<-opt$readthroughFilter
expressionFilter<-opt$expressionFilter
artifactFilter<-opt$artifactFilter
clinicalFile<-opt$clinicalFile
junctionReadCountFilter<-opt$junctionReadCountFilter
spanningFragCountFilter<-opt$spanningFragCountFilter
referenceFolder<-opt$referenceFolder
readingFrameFilter<-opt$readingFrameFilter

# standardFusioncallsSTARFusion<-readr::read_tsv(standardFusionFileSTARFusion)
# standardFusioncallsArriba<-readr::read_tsv(standardFusionFileArriba)
# 
# # combine callers
# standardFusioncalls<-rbind(standardFusioncallsArriba,standardFusioncallsSTARFusion)

standardFusioncalls<-lapply(standardFusionFiles,function(x) readr::read_tsv(x))
standardFusioncalls <- plyr::ldply(standardFusioncalls, data.frame) 


#formatting dataframe for filtering
standardFusioncalls<-standardFusioncalls %>%
  # to obtain geneA and geneB for gene search below
  bind_cols(colsplit(standardFusioncalls$FusionName, pattern = '--', names = c("GeneA","GeneB"))) %>%
  # Intergenic fusion will have Gene1A,Gene2A,Gene1B,Gene2B
  separate(GeneA,sep = "/",into =c("Gene1A","Gene2A"),remove=FALSE) %>%
  separate(GeneB,sep = "/",into =c("Gene1B","Gene2B"),remove=FALSE) %>%
  #remove distance to fusion breakpoint from gene names in intergenic fusion
  mutate(Gene1A=gsub("[(].*","",Gene1A),
         Gene2A=gsub("[(].*","",Gene2A),
         Gene1B=gsub("[(].*","",Gene1B),
         Gene2B=gsub("[(].*","",Gene2B)) %>%
  as.data.frame()


############# artifactFiltering #############


fusion_filtering_QC<-function(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,readthroughFilter=readthroughFilter){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param expressionMatrix A rds with expression data (RSEM/TPM counts) for the same cohort/GTEX
  # @param readingFramFilter A regex to capture readingframe (eg. in-frame|frameshift|other)
  # @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column: "annots"
  # @param junctionReadCountFilter An integer threshold for JunctionReadCount
  # @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
  # @param expressionFilter An integer threshold for TPM/FPKM filter
  # @return Standardized fusion calls filtered to pass QC and remove calls with insufficient read-support and annotation red-flags

  if( readthroughFilter & any(grepl('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE))){
    # Gather read throughs from standardized fusion calls
    rts <- standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),c("FusionName","annots")]
    if(length(grep("mitelman",rts$annots,ignore.case = TRUE))>0){
      #dont remove if fusion in mitelman (cancer fusion specific fusion database)
      rts <- rts[-grep("mitelman",rts$annots,ignore.case = TRUE),"FusionName"]
    }else{
      rts<-as.character(rts$FusionName)
    }
    # Reverse of read throughs to capture
    rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
    # Combine read through and reverse fusion genes
    rts <- unique(c(rts, rts.rev))
    # remove read throughs even if distance is not same in intergenic fusions
    rts<-unlist(lapply(rts,function(x) rm_between(x, "(", ")", extract = F)))
    rts<-data.frame("readThroughs"=rts)
    standardFusioncalls<-standardFusioncalls[-which(unlist(lapply(standardFusioncalls$FusionName,function(x) rm_between(x, "(", ")", extract = F))) %in% rts$readThroughs),]
  }

  if( !missing(readingFrameFilter ) ){
    # Error handling
    readingFrameTypes<-c("in-frame","frameshift","other")
    standardFusioncalls <- standardFusioncalls[grep(readingFrameFilter,standardFusioncalls$Fusion_Type),]
    if(!all(readingFrameTypes %in% standardFusioncalls$Fusion_Type)){
      warning(paste("No fusion calls with readingframe:",readingFrameTypes[-which(readingFrameTypes %in% standardFusioncalls$Fusion_Type)]))
    }
  }

  if( !missing(artifactFilter) & any(grepl(artifactFilter,standardFusioncalls$annots)) ) {
    # Error handling
    artifactFilterTypes<-unlist(strsplit(artifactFilter,"|",fixed=TRUE))
    if(any(unlist(lapply(artifactFilterTypes, function(x) !any(str_detect(as.character(standardFusioncalls$annots),x)))))){
      artifactFilterTypesNotFound<-artifactFilterTypes[unlist(lapply(artifactFilterTypes, function(x) !any(str_detect(as.character(standardFusioncalls$annots),x))))]
      warning(paste("No fusion calls with annotation:",artifactFilterTypesNotFound))
    }
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
QCFiltered<-fusion_filtering_QC(standardFusioncalls = standardFusioncalls,junctionReadCountFilter = junctionReadCountFilter,spanningFragCountFilter = spanningFragCountFilter,readingFrameFilter = readingFrameFilter,artifactFilter = artifactFilter,readthroughFilter = readthroughFilter)

saveRDS(QCFiltered,paste0(opt$outputFile,"_QC_filtered.RDS"))


############################################

########### Expression filtering ###########

# load expressionMatrix RDS for expression based filtering for less than given threshold
expressionMatrix<-readRDS(expressionMatrix)
# find the list of cohorts and sample type of interest
matched_samples <- read.delim(clinicalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(cohort %in% cohortInterest) %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")
# filter the expression to only the ones that are in the cohort and sample type of interest
expressionMatrix <- expressionMatrix %>%
  select(rownames(matched_samples)) %>%
  tibble::rownames_to_column("GeneSymbol")

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
  # Get the data into long format
  reshape2::melt(variable.name = "Sample",
                 value.name = "expression_value") %>%
  # Remove rows with expression that is too low
  dplyr::filter(expression_value > expressionFilter)

# Error handling
if (!all(fusion_sample_gene_df$Sample %in% expression_long_df$Sample)) {
  warning("Not all samples in expression file. Only returning fusions for samples in expressionMatrix.")
}

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
# column 1 as GeneName 2 source file 3 type; collapse to summarize type
geneListReferenceDataTab<-read.delim(file.path(referenceFolder,"genelistreference.txt"),stringsAsFactors = FALSE)
geneListReferenceDataTab<-geneListReferenceDataTab %>%
  # upper case because some genes have a/b/c etc
  mutate(Gene_Symbol=toupper(Gene_Symbol)) %>%
  dplyr::group_by(Gene_Symbol) %>%
  #collapse the gene type to have unique lines per gene
  dplyr::mutate(type = toString(type)) %>%
  dplyr::distinct(Gene_Symbol, type) %>% as.data.frame()

# column 1 as FusionName 2 source file 3 type; collapse to summarize type
fusionReferenceDataTab<-read.delim(file.path(referenceFolder,"fusionreference.txt"),stringsAsFactors = FALSE)
fusionReferenceDataTab<-fusionReferenceDataTab %>%
  dplyr::distinct(FusionName,type) %>% as.data.frame()


annotate_fusion_calls<-function(standardFusioncalls=standardFusioncalls,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab){
  # @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
  # @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
  # @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
  # @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder
  annotated_filtered_fusions<-standardFusioncalls %>%
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
  return(annotated_filtered_fusions)
}

# Annotate QC filtered fusion calls
annotated_filtered_fusions<-annotate_fusion_calls(standardFusioncalls = expression_filtered_fusions,geneListReferenceDataTab = geneListReferenceDataTab ,fusionReferenceDataTab = fusionReferenceDataTab )
saveRDS(annotated_filtered_fusions,paste0(opt$outputFile,"_QC_expression_filtered_annotated.RDS"))



############################################################################################

############################################################################################



