# K. S. Gaonkar 2019
# Annotates standardizes fusion calls from callers [STARfusion| Arriba] or QC filtered fusion
# calls with zscored expression value from either GTEx/cohort . The input should have the following standardized
# columns to run through this GTEx/cohort normalization function
# "Sample"  Unique SampleIDs used in your RNAseq dataset
# "FusionName" GeneA--GeneB ,or if fusion is intergenic then Gene1A/Gene2A--GeneB
#
# Optional
# "LeftBreakpoint" Genomic location of breakpoint on the left
# "RightBreakpoint" Genomic location of breakpoint on the right
# "Fusion_Type" Type of fusion prediction from caller "inframe","frameshift","other"
# "Caller" Name of caller used for fusion prediction
# "JunctionReadCount" synonymous with split reads, provide evidence for the specific breakpoint
# "SpanningFragCount" fragment alignment where either mate does not cross the junction
# "Confidence" Estimated confidence of call from caller, if not available value is "NA"
# "annots" Annotation from FusionAnnotator
#
# Command line arguments
# fusionfile	: Merged fusion calls from [STARfusion| Arriba] or QC filtered fusion
# zscoreFilter	: Zscore value to use as threshold for annotation of differential expression
# gtexMatrix : Rds will be saved as "normData" with expression data (FPKM) for samples (colnames) and
#                               GeneSymbol( rownames) to be used to normalize samples
# expressionMatrix : expression matrix (FPKM for samples that need to be zscore normalized)
# saveZscoredMatrix : path to save zscored matrix from GTEx normalization
# outputfile	: Output file is a standardized fusion call set with standard
# column headers with additional columns below from expression annotation
# zscore_Gene1A : zscore from GTEx/cohort normalization for Gene1A
# zscore_Gene1B : zscore from GTEx/cohort normalization for Gene1B
# zscore_Gene2A : zscore from GTEx/cohort normalization for Gene2A
# zscore_Gene2B : zscore from GTEx/cohort normalization for Gene2B
# note_expression_Gene1A : differentially expressed or no change in expression for Gene1A
# note_expression_Gene1B : differentially expressed or no change in expression for Gene1B
# note_expression_Gene2A : differentially expressed or no change in expression for Gene2A
# note_expression_Gene2B : differentially expressed or no change in expression for Gene2B




suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))

option_list <- list(
  make_option(c("-S", "--standardFusionCalls"),type="character",
              help="Standardized fusion calls (.RDS) "),
  make_option(c("-c", "--zscoreFilter"), type="integer",default=2,
              help="zscore value to use as threshold for annotation of differential expression"),
  make_option(c("-e","--expressionMatrix"),type="character",
              help="expression matrix (FPKM for samples that need to be zscore normalized .RDS)"),
  make_option(c("-s","--saveZscoredMatrix"),type="character",
              help="path to save zscored matrix from GTEx normalization (.RDS)"),
  make_option(c("-g","--gtexMatrix"),type="character",
              help="GTEx normalization FPKM expression data to compare (.rds)"),
  make_option(c("-o","--outputfile"),type="character",
              help="Standardized fusion calls with additional annotation from GTEx zscore comparison (.TSV)")
)

opt <- parse_args(OptionParser(option_list=option_list))
standardFusionCalls<-opt$standardFusionCalls
expressionMatrix<-opt$expressionMatrix
zscoreFilter<- opt$zscoreFilter
saveZscoredMatrix<-opt$saveZscoredMatrix
gtexMatrix<-opt$gtexMatrix


ZscoredAnnotation<-function(standardFusionCalls=standardFusionCalls,zscoreFilter=zscoreFilter,saveZscoredMatrix=saveZscoredMatrix,normData=normData,expressionMatrix=expressionMatrix){
  #  @param standardFusionCalls : Annotates standardizes fusion calls from callers [STARfusion| Arriba] or QC filtered fusion
  #  @param zscoreFilter : Zscore value to use as threshold for annotation of differential expression
  #  @results : expression_annotated_fusions is a standardized fusion call set with standard
              # column headers with additional columns below from expression annotation
              # zscore_Gene1A : zscore from GTEx normalization for Gene1A
              # zscore_Gene1B : zscore from GTEx normalization for Gene1B
              # zscore_Gene2A : zscore from GTEx normalization for Gene2A
              # zscore_Gene2B : zscore from GTEx normalization for Gene2B
              # note_expression_Gene1A : differentially expressed or no change in expression for Gene1A
              # note_expression_Gene1B : differentially expressed or no change in expression for Gene1B
              # note_expression_Gene2A : differentially expressed or no change in expression for Gene2A
              # note_expression_Gene2B : differentially expressed or no change in expression for Gene2B

  # expressionMatrix collapsed at gene level like Gtex max rowMean
  expressionMatrixMatched <- expressionMatrix %>%
    unique() %>%
    # means for each row per each gene_id
    dplyr::mutate(means = rowMeans(select(.,-GeneSymbol,-gene_id,-EnsembleID))) %>%
    # arrange descending mean
    arrange(desc(means)) %>%
    # to keep only first occurence ie. max rowMean per GeneSymbol
    distinct(GeneSymbol, .keep_all = TRUE) %>%
    ungroup() %>%
    dplyr::select(-means,-gene_id,-EnsembleID) %>%
    dplyr::filter(GeneSymbol %in% normData$GeneSymbol) %>%
    tibble::column_to_rownames("GeneSymbol")
  expressionMatrixMatched<-log2(expressionMatrixMatched+1)



  # gene matched
  # get log transformed GTEx matrix
  normData <- normData %>%
    tibble::column_to_rownames("GeneSymbol") %>%
    as.matrix()
  normData <- normData[rownames(expressionMatrixMatched), ]
  normData <- log2(normData + 1)

  #normData mean and sd
  normData_means <- rowMeans(normData, na.rm = TRUE)
  normData_sd <- apply(normData, 1, sd, na.rm = TRUE)
  # subtract mean
  expressionMatrixzscored <- sweep(expressionMatrixMatched, 1, normData_means, FUN = "-")
  # divide by SD
  expressionMatrixzscored <- sweep(expressionMatrixzscored, 1,normData_sd, FUN = "/") %>% mutate(GeneSymbol=rownames(.)) %>% na.omit()


  # To save GTEx/cohort scored matrix
  if(!missing(saveZscoredMatrix)){
    saveRDS(expressionMatrixzscored,saveZscoredMatrix)
  }

  # get long format to compare to expression
  expression_long_df <- expressionMatrixzscored %>%
    # Get the data into long format
    reshape2::melt(variable.name = "Sample",
                   value.name = "zscore_value")

    # fusion calls
  fusion_sample_gene_df <- standardFusionCalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(Sample, FusionName, Gene1A, Gene1B, Gene2A, Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                  key = gene_position, value = GeneSymbol) %>%
    # Remove columns without gene symbols
    dplyr::filter(GeneSymbol != "") %>%
    # This is for illustrations sake, only
    dplyr::arrange(Sample, FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct()

  expression_annotated_fusions <- fusion_sample_gene_df %>%
    # join the filtered expression values to the data frame keeping track of symbols
    # for each sample-fusion name pair
    dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol")) %>%
    # for each sample-fusion name pair, are all genes under the expression threshold?
    dplyr::group_by(FusionName, Sample) %>%
    dplyr::select(FusionName, Sample,zscore_value,gene_position) %>%
    # cast to keep zscore and gene position
    reshape2::dcast(FusionName+Sample ~gene_position,value.var = "zscore_value") %>%
    # get annotation from z score
    dplyr::mutate(note_expression_Gene1A = ifelse((Gene1A>zscoreFilter| Gene1A< -zscoreFilter),"differentially expressed","no change"),
                  note_expression_Gene1B = ifelse((Gene1B>zscoreFilter| Gene1B< -zscoreFilter),"differentially expressed","no change"),
                  note_expression_Gene2A = ifelse((Gene2A>zscoreFilter| Gene2A< -zscoreFilter),"differentially expressed","no change"),
                  note_expression_Gene2B = ifelse((Gene2B>zscoreFilter| Gene2B< -zscoreFilter),"differentially expressed","no change")) %>%
    dplyr::rename(zscore_Gene1A=Gene1A,
                  zscore_Gene1B=Gene1B,
                  zscore_Gene2A=Gene2A,
                  zscore_Gene2B=Gene2B) %>%
    # unique FusionName-Sample rows
    # use this to filter the QC filtered fusion data frame
    dplyr::inner_join(standardFusionCalls, by = c("FusionName", "Sample")) %>%
    dplyr::distinct()

  return (expression_annotated_fusions)
}

# example run using GTEx
# load GTEx data
normData<-read_tsv(gtexMatrix) %>% rename(GeneSymbol=gene_symbol)

# load standardaized fusion calls cohort
standardFusionCalls<-readRDS(standardFusionCalls)

# load expression Matrix for cohort
expressionMatrix<-readRDS(expressionMatrix)
expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))

# for cohort level run the expressionMatrix divided by broad_histology will be provided to the function as normData and expressionMatrix

GTExZscoredAnnotation_filtered_fusions<- ZscoredAnnotation(standardFusionCalls ,zscoreFilter,normData=normData,expressionMatrix = expressionMatrix)
write.table(GTExZscoredAnnotation_filtered_fusions,paste0(opt$outputfile,"_GTExComparison_annotated.RDS"))





