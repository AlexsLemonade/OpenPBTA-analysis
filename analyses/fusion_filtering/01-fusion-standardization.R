# K. S. Gaonkar 2019
# Standardizes fusion calls from callers [STARfusion| Arriba]. The output will
# have the following columns
# "Sample"  Unique SampleIDs used in your RNAseq dataset
# "LeftBreakpoint" Genomic location of breakpoint on the left
# "RightBreakpoint" Genomic location of breakpoint on the right
# "FusionName" GeneA--GeneB ,or if fusion is intergenic then Gene1A/Gene2A--GeneB
# "Fusion_Type" Type of fusion prediction from caller "inframe","frameshift","other"
# "Caller" Name of caller used for fusion prediction
# "JunctionReadCount" synonymous with split reads, provide evidence for the specific breakpoint
# "SpanningFragCount" fragment alignment where either mate does not cross the junction
# "Confidence" Estimated confidence of call from caller, if not available value is "NA"
# "annots" Annotation from FusionAnnotator
#
# Command line arguments
# fusionfile	: Merged fusion calls from [STARfusion| Arriba]
# For STARfusion default columns must be provided .For Arriba along with default
# columns additional columns from Fusion Annotator [fusion_name= GeneA--GeneB;
# annots=Annotation from FusionAnnotator] should be appended to the right of
# original arriba data frame
#
# caller	: Caller used for the input merged fusion calls [STARfusion| Arriba]
#
# outputfile	: Output file is a standardized fusion call set with standard
# column headers used in the filtering scripts. Column order:
# 1:  'LeftBreakpoint
# 2:  'RightBreakpoint'
# 3:  'FusionName'
# 4:  'Sample',
# 5:  'Caller'
# 6:  'Fusion_Type'
# 7:  'JunctionReadCount'
# 8:  'SpanningFragCount'
# 9:  'Confidence'
# 10: 'annots'



suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
     make_option(c("-f", "--fusionfile"),type="character",
        help="Merged fusion calls from [STARfusion | Arriba]"),
     make_option(c("-c", "--caller"), type="character",
        help="Caller type [STARfusion|| Arriba]"),
     make_option(c("-o","--outputfile"),type="character",
	help="Standardized fusion calls from [STARfusion | Arriba] (.TSV)")
    )

# Get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))
inputfile <- opt$fusionfile
caller <- opt$caller
# Converting caller ID to call Uppercase for error handling
caller<-toupper(caller)
outputfile <- opt$outputfile

fusion_calls <- read.delim(inputfile,stringsAsFactors=FALSE)
# To have a general column with unique IDs associated with each sample
fusion_calls$Sample <- fusion_calls$tumor_id
fusion_calls$Caller <- caller



standard_fusion <- function(fusion_calls=fusion_calls,caller=caller) {
# @param fusion_calls A dataframe from star fusion or arriba (more callers to be added)
# @param caller string options STARfusion/arriba
# @return Standardized fusion calls ready for filtering

  if( caller == "STARFUSION"){
    fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$LeftBreakpoint)
    fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$RightBreakpoint)
    # Standardizing fusion type annotation
    fusion_calls <- fusion_calls %>%
    rename(Fusion_Type = PROT_FUSION_TYPE) %>%
    dplyr::mutate(Fusion_Type = case_when(
    Fusion_Type == "INFRAME" ~ "in-frame",
    Fusion_Type == "FRAMESHIFT" ~ "frameshift",
    TRUE ~ "other"
    ))
    # Adding Confidence column
    fusion_calls$Confidence<-"NA"
  }
  else if( caller == "ARRIBA"){
    fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$breakpoint1)
    fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$breakpoint2)
    # Standardizing fusion type annotation
    fusion_calls <- fusion_calls %>%
    rename(Fusion_Type = reading_frame) %>%
    dplyr::mutate(Fusion_Type = case_when(
    !Fusion_Type %in% c("out-of-frame", "in-frame") ~ "other",
    Fusion_Type == "out-of-frame" ~ "frameshift",
    TRUE ~ "in-frame"
    ))
    colnames(fusion_calls)[which(colnames(fusion_calls)=="X..")] <- "annots"
    # Intergenic gene fusion breakpoints in arriba are annotated as
    # "gene1A,gene1B". As comma is used as a common delimiter in files changing
    # it to "/"
    fusion_calls$FusionName <- paste0(gsub("," , "/" , fusion_calls$gene1) ,"--" ,
    gsub("," , "/" , fusion_calls$gene2))
    # JunctionReadCount is equivalent to split reads in Arriba. Arriba however
    # provides split_reads1 and split_reads2 to provide information of reads
    # anchoring in gene1 or gene2
    fusion_calls$JunctionReadCount <- fusion_calls$split_reads1+fusion_calls$split_reads2
    # SpanningFragCount is equivalent to discordant_mates in Arriba
    fusion_calls$SpanningFragCount <- fusion_calls$discordant_mates
    fusion_calls$Confidence <- fusion_calls$confidence
    # To rename Arriba dataframe's last column which should be output from FusionAnnotator
  } else {
    stop(paste(caller, "is not a supported caller string."))
  }

  #Get standard columns for filtering
	standard_calls <- unique(fusion_calls[,c('LeftBreakpoint','RightBreakpoint','FusionName' , 'Sample' , 'Caller' ,
	'Fusion_Type' , 'JunctionReadCount' , 'SpanningFragCount' , 'Confidence' ,
	'annots')])
  return(standard_calls)


}

# Standardized fusion calls ready for filtering
output <- standard_fusion(fusion_calls = fusion_calls , caller = caller)

saveRDS(output , outputfile)
write.table(output , file=outputfile , sep="\t" , row.names = FALSE , col.names = TRUE,quote=FALSE)
