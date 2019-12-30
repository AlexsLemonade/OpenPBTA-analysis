# Author: Krutika Gaonkar
#
# Read in concensus snv calls to gather alterations in TP53 and NF1 
# to evaluate classifier
# @params snvConcensus multi-caller concensus snv calls
# @params clincalFile clinical file: pbta-histologies.tsv 
# @params outputFolder output folder for alteration file


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))


option_list <- list(
  make_option(c("-s", "--snvConsensus"),type="character",
              help="Consensus snv calls (.tsv) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-o","--outputFolder"),type="character",
              help="output folder for results ")
)

opt <- parse_args(OptionParser(option_list=option_list))
snvConsensusFile<-opt$snvConsensus
clinicalFile<-opt$clinicalFile
outputFolder<-opt$outputFolder

# read in snv
consensus_snv<-read_tsv(snvConsensusFile)
TP53_alt<-consensus_snv %>% filter(Hugo_Symbol=="TP53")
NF1_alt<-consensus_snv %>% filter(Hugo_Symbol=="NF1")

# clinical file
clinical<-read_tsv(clinicalFile)
clinical<-clinical %>% filter(!sample_type=="Normal" & !experimental_strategy == "RNA-Seq")
# fix PNOC sample ids?
clinical[grep("-T",clinical$sample_id),"sample_id"]<-lapply(as.vector(clinical[grep("-T",clinical$sample_id),"sample_id"]),function(x) gsub("-T$","",gsub("^.*-T-","",gsub("[.].*","",x))))


# remove non-coding mutations
nonCodingVariant<-c("Intron","3'Flank","5'Flank","3'UTR","5'UTR","Silent")
TP53_alt<-TP53_alt %>% filter(!Variant_Classification %in% nonCodingVariant) %>% full_join(clinical,by=c("Tumor_Sample_Barcode"="Kids_First_Biospecimen_ID")) 
TP53_alt$Hugo_Symbol[is.na(TP53_alt$Hugo_Symbol)]="No_TP53_alt"

# remove non-coding mutations
NF1_alt<-NF1_alt %>% filter(!Variant_Classification %in% nonCodingVariant) %>% full_join(clinical,by=c("Tumor_Sample_Barcode"="Kids_First_Biospecimen_ID")) 
NF1_alt$Hugo_Symbol[is.na(NF1_alt$Hugo_Symbol)]="No_NF1_alt"


# save TP53 and NF1 SNV alterations
TP53_NF1_alt<-rbind(TP53_alt,NF1_alt) 
write.table(TP53_NF1_alt,file.path(outputFolder,"TP53_NF1_snv_alteration.tsv"),sep="\t",quote=FALSE,row.names = FALSE)





