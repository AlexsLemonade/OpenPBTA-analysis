# K. S. Gaonkar 2019
# Standardizes fusion calls from callers [STARfusion|| Arriba]. The output will have the following columns
# "Sample"  Unique SampleIDs used in your RNAseq dataset
# "FusionName" GeneA--GeneB if fusion is intergenic then Gene1A/Gene2A--GeneB  	
# "Fusion_Type" Type of fusion prediction from caller "inframe","frameshift","other"
# "Caller" Name of caller used for fusion prediction
# "JunctionReadCount" synonymous with split reads, provide evidence for the specific breakpoint 
# "SpanningFragCount" fragment alignment where either mate does not cross the junction
# "Confidence" Estimated confidence of call from caller, if not available value is "NA"
# "annots" Annotation from FusionAnnotator
# 
# Command line arguments
# fusionfile : 	Merged fusion calls from [STARfusion|| Arriba]
#		For STARfusion default columns must be provided
#		For Arriba along with default columns additional 
#		columns from Fusion Annotator [fusion_name= GeneA--GeneB; 
#		annots=Annotation from FusionAnnotator] are appended to the 
#		right of original arriba data frame 
#
# caller :	Caller used for the input merged fusion calls [STARfusion|| Arriba]
#

library("optparse")

option_list <- list( 
     make_option(c("-f", "--fusionfile"),type="character",
        help="Merged fusion calls from [STARfusion || Arriba]"),
     make_option(c("-c", "--caller"), type="character", 
        help="Caller type [STARfusion|| Arriba]"),
     make_option(c("-o","--outputfile"),type="character",
	help="Standardized fusion calls from [STARfusion || Arriba] (.RDS)")	
    )

# get command line options, if help option encountered print help and exit,

opt <- parse_args(OptionParser(option_list=option_list))

inputfile<-opt$fusionfile
caller<-opt$caller
outputfile<-opt$outputfile

fusion_calls<-read.delim(opt$fusionfile,stringsAsFactors=FALSE)
#To have a general column with unique IDs associated with each sample
fusion_calls$Sample<-fusion_calls$tumor_id



standard_fusion <- function(fusion_calls=fusion_calls,caller=caller) {
# @param fusion_calls A dataframe from star fusion or arriba (more callers to be added)
# @param caller string options STARfusion/arriba
# @return Standardized fusion calls ready for filtering
  
  if(caller=="STARfusion"){
  	fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$LeftBreakpoint)
  	fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$RightBreakpoint)
  	fusion_calls$Fusion_Type<-fusion_calls$PROT_FUSION_TYPE
  	fusion_calls$Fusion_Type<-gsub("INFRAME","in-frame",fusion_calls$Fusion_Type)
  	fusion_calls$Fusion_Type<-gsub("FRAMESHIFT","frameshift",fusion_calls$Fusion_Type)
  	fusion_calls$Fusion_Type[-which(fusion_calls$Fusion_Type %in% c("in-frame","frameshift"))] <- 'other'
  	fusion_calls$Caller <- 'STARFusion'
  	fusion_calls$FusionName<-fusion_calls$FusionName
  	fusion_calls$Confidence<-"NA"
  	standard_calls <- unique(fusion_calls[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence','annots')])
  	return(standard_calls)

  }
  if( caller=="Arriba"){
	fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$breakpoint1)
    	fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$breakpoint2)
    	fusion_calls$Fusion_Type<-fusion_calls$reading_frame
    	fusion_calls$Fusion_Type<-gsub("out-of-frame","frameshift",fusion_calls$Fusion_Type)
    	fusion_calls$Fusion_Type[-which(fusion_calls$Fusion_Type %in% c("in-frame","frameshift"))] <- 'other'
    	fusion_calls$Caller <- 'arriba'
    	fusion_calls$FusionName <-paste0(gsub(",","/",fusion_calls$gene1),"--",gsub(",","/",fusion_calls$gene2))
    	fusion_calls$JunctionReadCount<-fusion_calls$split_reads1+fusion_calls$split_reads2
    	fusion_calls$SpanningFragCount<-fusion_calls$discordant_mates
    	fusion_calls$Confidence<-fusion_calls$confidence
	colnames(fusion_calls)[27]<-"annots"
	standard_calls <- unique(fusion_calls[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence','annots')])
    	return(standard_calls)
  }

}

output<-standard_fusion(fusion_calls=fusion_calls,caller=caller)

saveRDS(output,outputfile)



