# Author: Krutika Gaonkar
#
# Read in concensus snv calls to gather alterations in TP53 and NF1 
# to evaluate classifier
# @params snvConcensus multi-caller concensus snv calls
# @params clincalFile clinical file: pbta-histologies.tsv 
# @params outputFolder output folder for alteration file
# @params gencode cds bed file from gencode 


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("GenomicRanges"))


option_list <- list(
  make_option(c("-s", "--snvConsensus"),type="character",
              help="Consensus snv calls (.tsv) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-o","--outputFolder"),type="character",
              help="output folder for results "),
  make_option(c("-g","--gencode"),type="character",
              help="cds gencode bed file")
)

opt <- parse_args(OptionParser(option_list=option_list))
snvConsensusFile<-opt$snvConsensus
clinicalFile<-opt$clinicalFile
outputFolder<-opt$outputFolder
gencodeBed<-opt$gencode

# read in snv
consensus_snv<-read_tsv(snvConsensusFile)
# gencode cds region
gencode<-read.delim(gencodeBed,stringsAsFactors = F,header = F)
gencode_gr <- with(gencode, GRanges(V1, IRanges(V2, V3)))

maf_to_granges<-function(maf=maf){
  gr<-with(maf,GRanges(maf$Chromosome, IRanges(maf$Start_Position, maf$End_Position), Variant_Classification=Variant_Classification, Tumor_Sample_Barcode=Tumor_Sample_Barcode, Hugo_Symbol=Hugo_Symbol))
  return(gr)
}

# TP53 alterations 
TP53_alt<-consensus_snv %>% filter(Hugo_Symbol=="TP53")
TP53_gr<-maf_to_granges(TP53_alt)

#NF1 alterations
NF1_alt<-consensus_snv %>% filter(Hugo_Symbol=="NF1")
NF1_gr<-maf_to_granges(NF1_alt)

intersect_bed <- function(a, b){
  my_hit <- findOverlaps(a, b)
  my_df  <- cbind(as.data.frame(a[queryHits(my_hit)]),
                  as.data.frame(b[subjectHits(my_hit)]))
  return (my_df)
}

TP53_cds_overlap<-intersect_bed(TP53_gr,gencode_gr)
NF1_cds_overlap<-intersect_bed(NF1_gr,gencode_gr)


# clinical file
clinical<-read_tsv(clinicalFile)
clinical<-clinical %>% filter(!sample_type=="Normal" & !experimental_strategy == "RNA-Seq")
# fix PNOC sample ids?
clinical[grep("-T",clinical$sample_id),"sample_id"]<-lapply(as.vector(clinical[grep("-T",clinical$sample_id),"sample_id"]),function(x) gsub("-T$","",gsub("^.*-T-","",gsub("[.].*","",x))))


# remove non-coding mutations
nonCodingVariant<-c("Intron","3'Flank","5'Flank","3'UTR","5'UTR","Silent")

# merge TP53 and NF1 mutations (remove subject seqname start end)
TP53_NF1_cds_alt<-rbind(TP53_cds_overlap[,-c(9:13)],NF1_cds_overlap[,-c(9:13)]) 

TP53_NF1_cds_alt<-TP53_NF1_cds_alt  %>% filter(!Variant_Classification %in% nonCodingVariant) %>% full_join(clinical,by=c("Tumor_Sample_Barcode"="Kids_First_Biospecimen_ID")) 
TP53_NF1_cds_alt$Hugo_Symbol[is.na(TP53_NF1_cds_alt$Hugo_Symbol)]="No_TP53_NF1_alt"


# save TP53 and NF1 SNV alterations
write.table(TP53_NF1_cds_alt,file.path(outputFolder,"TP53_NF1_snv_alteration.tsv"),sep="\t",quote=FALSE,row.names = FALSE)





