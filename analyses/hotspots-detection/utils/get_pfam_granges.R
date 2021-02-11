#' Given pfam id(s) and an ensembl mart 
#' creates a grange for the pfam domains
#' @param pfam_id pfam id(s) to filter 
#' ensembl 
#' @param ensembl a S4 class from biomaRt::useMart() call
#' 
#'@return granges for genomic locations of given pfam id(s)
#'
#'

get_pfam_granges <- function(pfam_id,ensembl){
  # hsapiens_gene_ensembl == Human genes (GRCh38.p13) 10th Feb 2021 
  
  pfamDataBioMart<-getBM(attributes = c('pfam','hgnc_symbol'),  
                         filters = "pfam",
                         values = pfam_id,
                         mart = ensembl)
  # Pfam file downloaded from
  # contains pfam IDs and associated names
  pfam<-read_tsv(file.path("input","pfamDesc.txt.gz"),col_names = FALSE) %>% as.data.frame()
  colnames(pfam) <- c("ID", "NAME", "DESC")
  
  # UCSC file downloaded from
  # contains pfam names and chromosome positions
  pfam_location <- read.delim(file.path("input","ucscGenePfam.txt.gz"), header=F)
  colnames(pfam_location) <- c("BIN", "CHROM", "CHROM_START", "CHROM_END", "NAME", "SCORE",
                               "STRAND", "THICK_START", "THICK_END", "RESERVED", "BLOCK_COUNT",
                               "BLOCK_SIZES", "CHROM_STARTS")
  
  # merge using pfam names tp get pfam ID
  pfamDescLoc <- pfam %>%
    left_join(pfam_location, by="NAME") %>%
    dplyr::select("ID","NAME","DESC", "CHROM", "CHROM_START", "CHROM_END", "STRAND")
  
  # get pfam description and name for pfam id
  pfamGene <- pfamDataBioMart %>% left_join(pfamDescLoc, by=c("pfam"="ID"))
  
  # granges
  pfamGene_gr <- makeGRangesFromDataFrame(df = pfamGene,keep.extra.columns = TRUE,seqnames.field = "CHROM",start.field = "CHROM_START",end.field = "CHROM_END",strand.field = "STRAND")
}
