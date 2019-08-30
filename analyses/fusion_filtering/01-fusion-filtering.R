library("reshape2")
library("data.table")
library("plyr")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(R.utils)

# star fusion
sf<-read.delim(gzfile(paste0(file.path("../../data/pbta-fusion-starfusion.tsv.gz"))),stringsAsFactors=F,,header=T,sep="\t")
#head(sf)
sf<-sf[-which(sf$SpanningFragCount-sf$JunctionReadCount >10|sf$JunctionReadCount==0|sf$LargeAnchorSupport == "NO_LDAS"),]
sf$LeftBreakpoint <- gsub('^chr','',sf$LeftBreakpoint)
sf$RightBreakpoint <- gsub('^chr','',sf$RightBreakpoint)
sf$Fusion_Type<-sf$PROT_FUSION_TYPE
sf$Fusion_Type[which(sf$PROT_FUSION_TYPE=="INFRAME")] <- 'in-frame'
sf$Fusion_Type[grep("FRAMESHIFT",sf$PROT_FUSION_TYPE)] <- 'frameshift'
sf$Fusion_Type[-which(sf$PROT_FUSION_TYPE %in% c("INFRAME","FRAMESHIFT"))] <- 'other'
sf$Caller <- 'STARFusion'
sf$Sample <-sf$tumor_id
sf$FusionName<-sf$FusionName
sf$Confidence<-"NA"
head(sf)
sf.total <- unique(sf[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])

#sf.total<-aggregate(list(sf.total$JunctionReadCount,sf.total$JunctionReadCount), list(sf.total$Sample,sf.total$FusionName,sf.total$Caller,sf.total$Fusion_Type, sf.total$Confidence,sf.total$SpanningFragCount) ,sum)
#colnames(sf.total)<-c('Sample','FusionName','Caller','Fusion_Type','Confidence','SpanningFragCount','JunctionReadCount')



#IGLC6 is expressed in rnaseq
sf.total$FusionName<-sub("IGL@","IGLC6",sf.total$FusionName)
sf.total$FusionName<-sub("IGL-@","IGLC6",sf.total$FusionName)

# arriba fusion
ar <- read.delim(gzfile(paste0(file.path("../../data/pbta-fusion-arriba.tsv.gz"))),stringsAsFactors=F,,header=T,sep="\t")
#remove false positives through supporting reads
ar<-ar[-(which(ar$discordant_mates-(ar$split_reads1+ar$split_reads2) >10 |ar$split_reads1+ar$split_reads2==0)),]
ar$LeftBreakpoint <- gsub('^chr','',ar$breakpoint1)
ar$RightBreakpoint <- gsub('^chr','',ar$breakpoint2)
ar$Fusion_Type<-ar$reading_frame
ar$Fusion_Type[grep("out-of-frame",ar$Fusion_Type)]<-"frameshift"
ar$Fusion_Type[-which(ar$Fusion_Type %in% c("in-frame","frameshift"))] <- 'other'

ar$Caller <- 'arriba'
ar$Sample <-ar$tumor_id
ar$FusionName <-paste0(gsub(",","/",ar$gene1),"--",gsub(",","/",ar$gene2))
ar$SpanningFragCount<-ar$discordant_mates
ar$JunctionReadCount<-ar$split_reads1+ar$split_reads2
ar$Confidence<-ar$confidence
ar.total <- unique(ar[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])

#ar.total<-aggregate(list(ar.total$JunctionReadCount,ar.total$SpanningFragCount), list(ar.total$Sample,ar.total$FusionName,ar.total$Caller,ar.total$Fusion_Type, ar.total$Confidence) ,sum)
#colnames(ar.total)<-c('Sample','FusionName','Caller','Fusion_Type','Confidence','SpanningFragCount','JunctionReadCount')



#reassign IGL@ gene name
ar.total$FusionName<-sub("IGL@","IGLC6",ar.total$FusionName)
ar.total$FusionName<-sub("IGL-@","IGLC6",ar.total$FusionName)


#gather read throughs

sf.rt <- unique(sf[grep('readthrough|neighbors|GTEx_Recurrent|BodyMap|DGD_PARALOGS|HGNC_GENEFAM|ConjoinG',sf$annots,ignore.case=TRUE),])
#https://arriba.readthedocs.io/en/latest/interpretation-of-results/#frequent-types-of-false-positives
#not filtered for pcr_fusions
#only removing read-throughs 
ar.rt <- unique(ar[grep("read-through",ar$type),])
ar.rt <- rbind(ar.rt,unique(ar[grep("readthrough|Normal|neighbors|GTEx_Recurrent|BodyMap|DGD_PARALOGS|HGNC_GENEFAM|ConjoinG",ar[,27],ignore.case=TRUE),]))
rts <- unique(c(sf.rt$FusionName, ar.rt$FusionName))
rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
rts <- unique(c(rts, rts.rev))


#merge callers
all.callers<-rbind(ar.total, sf.total)

#histology
clin<-read.delim(paste0(file.path("../../data/pbta-histologies.tsv")),stringsAsFactors = F,sep="\t")
clin$sample_id<-clin$Kids_First_Biospecimen_ID

#merge callers and clinical information
all.callers<-merge(all.callers, clin, by.x="Sample", by.y="sample_id")
all.callers <- unique(all.callers[,c('Sample','FusionName','Caller','Fusion_Type','broad_histology','JunctionReadCount','SpanningFragCount','Confidence')])
colnames(all.callers)<-c("Sample","Fused_Genes","Caller","Fusion_Type","Histology.Broad","JunctionReadCount","SpanningFragCount","Confidence")


# remove read-throughs
all.callers <- all.callers[-which(all.callers$Fused_Genes %in% rts),]
nrow(all.callers)


#public databases to gather putative driver fusion calls
tsgs<-read.delim('../../scratch/fusion_filtering_pipeline/references/tsgs.txt', header = T,stringsAsFactors = F)
onco<-read.delim('../../scratch/fusion_filtering_pipeline/references/allOnco_Feb2017.tsv', header = T,stringsAsFactors = F)
tcga<-read.delim('../../scratch/fusion_filtering_pipeline/references/pancanfus.txt', header = T,stringsAsFactors = F)
cosmic<-read.delim('../../scratch/fusion_filtering_pipeline/references/Cosmic_gene_census.csv',header=T,stringsAsFactors=F,sep=",")
tcga$TCGA_fusions<-paste(tcga$Gene_A,tcga$Gene_B,sep="--")
head(tcga)

#capture tsgs,onco,tcga fusion genes
genes.to.search <- c(paste0('^',tsgs$GeneSymbol,'-'), paste0('-',tsgs$GeneSymbol,'$'),paste0('^',onco$symbol,'-'), paste0('-',onco$symbol,'$'),paste0('^',cosmic$Gene.Symbol,'-'),paste0('-',cosmic$Gene.Symbol,'$'))
fusion.to.search<-tcga$TCGA_fusions
tsgs_onco_fusion<-all.callers[unlist(lapply(genes.to.search,function(x) grep(x, all.callers$Fused_Genes))),]
tsgs_onco_fusion<-rbind(tsgs_onco_fusion,all.callers[unlist(lapply(fusion.to.search,function(x) grep(x, all.callers$Fused_Genes))),])
tsgs_onco_fusion$note<-"Found in onco/tsgs/cosmic/tcga list"
to.add<-tsgs_onco_fusion
to.add$Caller_type<-paste(to.add$Caller,to.add$Fusion_Type,to.add$Confidence,sep="_")



# Gene fusion should be in-frame
# Called by at least 2 callers
all.callers.summary <- all.callers %>% 
  filter(Fusion_Type != "other") %>%
  group_by(Fused_Genes, Sample, Histology.Broad ) %>% 
  unique() %>%
  mutate(Caller = toString(Caller), caller.count = n()) %>%
  filter(caller.count >= 2) %>% 
  select(-caller.count, -Caller, -Fusion_Type) %>%
  unique() %>%
  as.data.frame()

print("caller count")
all.callers.summary$note<-"Called by both callers"
nrow(all.callers.summary)

# or found in at least 2 samples of the same histology 
sample.count <- all.callers %>% 
  filter(Fusion_Type != "other") %>%
  group_by(Fused_Genes, Histology.Broad) %>% 
  unique() %>%
  mutate(sample.count = n()) %>%
  filter(sample.count > 1) %>%
  select(-Caller, -sample.count, -Fusion_Type) %>%
  unique() %>%
  as.data.frame() 
length(unique(sample.count$Fused_Genes))
sample.count$note<-"found in atleast 2 samples in same histology"

# or GeneB or GeneA gene recurrently fused within a histology (>= 5 genes)
rec <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
rec2 <- rec %>% group_by(Histology.Broad) %>% 
  select(Histology.Broad,GeneA,GeneB) %>% 
  unique() %>% group_by(Histology.Broad, GeneA) %>% 
  summarise(GeneA.ct = n()) %>%
  filter(GeneA.ct >= 5) %>% as.data.frame()
rec3 <- rec %>% group_by(Histology.Broad) %>% 
  select(Histology.Broad,GeneA,GeneB) %>% 
  unique() %>% group_by(Histology.Broad, GeneB) %>% 
  summarise(GeneB.ct = n()) %>%
  filter(GeneB.ct >= 5) %>% as.data.frame()
rec2 <- merge(rec2, rec, by = c('GeneA','Histology.Broad'))
rec3 <- merge(rec3, rec, by = c('GeneB','Histology.Broad'))
rec2 <- unique(rec2[,c("Sample","Fused_Genes","Histology.Broad","JunctionReadCount","SpanningFragCount","Confidence")])
rec3 <- unique(rec3[,c("Sample","Fused_Genes","Histology.Broad","JunctionReadCount","SpanningFragCount","Confidence")])
res <- unique(rbind(rec2, rec3))

res$note<-"recurrently fused in a histology"


# merge these 
total <- unique(rbind(all.callers.summary, sample.count, res))

#add TF and Kinase gene fusions
curatedtf <- read.delim(paste0('../../scratch/fusion_filtering_pipeline/','references/curatedTF_attribute_list_entries.txt'), header = T)
predictedtf<-read.delim(paste0('../../scratch/fusion_filtering_pipeline/','references/predictedTF_attribute_list_entries.txt'), header = T)
tf<-rbind(curatedtf,predictedtf)


kinase <- read.delim(paste0('../../scratch/fusion_filtering_pipeline/','/references/Kincat_Hsap.08.02.txt'),sep="\t",stringsAsFactors = F)
kinase <- kinase[-which(kinase$Entrez_Symbol == ""),]

genes.to.search <- c(paste0('^',kinase$Entrez_Symbol,'-'), paste0('-',kinase$Entrez_Symbol,'$'),paste0('^',tf$GeneSym,'-'), paste0('-',tf$GeneSym,'$'))
tf_kinase<-all.callers[unlist(lapply(genes.to.search,function(x) grep(x, all.callers$Fused_Genes))),]
tf_kinase$note<-"Found in TF/Kinase list"

colnames(tf_kinase)

total<-rbind(total,tf_kinase[,c("Sample","Fused_Genes","Histology.Broad","JunctionReadCount","SpanningFragCount","Confidence","note")])



# remove GeneA == GeneB fusion calls
#total <- cbind(total, colsplit(total$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
#total<- total[-which(total$GeneA == total$GeneB),]
#total<-total[,c("Fused_Genes","Sample","Histology.Broad","note","JunctionReadCount","SpanningFragCount")]


print("total number of calls")
head(total)

# remove fusions that are in > 1 histology 
hist.count <- total %>% 
  select(Fused_Genes, Histology.Broad,JunctionReadCount,SpanningFragCount) %>%
  unique() %>%
  group_by(Fused_Genes) %>%
  summarise(hist.count = n()) %>%
  filter(hist.count == 1)
total <- total[which(total$Fused_Genes %in% hist.count$Fused_Genes),]
length(unique(total$Fused_Genes))


#merge prioritised calls
total <- merge(total, all.callers, by = c('Sample','Fused_Genes','Histology.Broad','JunctionReadCount','SpanningFragCount','Confidence'))
total$Caller_type<-paste(total$Caller,total$Fusion_Type,total$Confidence,sep="_")


#filter for multi_fused
rec <- cbind(colsplit(total$Fused_Genes, pattern = '--', names = c("Gene1","Gene2")), total)
rec <- cbind(colsplit(rec$Gene2, pattern = '/', names = c("Gene2A","Gene2B")), rec)
rec <- cbind(colsplit(rec$Gene1, pattern = '/', names = c("Gene1A","Gene1B")), rec)

total<-rec
rec<-unique(rec[,c("Sample","Gene1","Gene1A","Gene1B","Gene2A","Gene2B","Gene2")])

rec[grep("\\(", rec$Gene1B),"Gene1B"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene1B),"Gene1B"])
rec[grep("\\(", rec$Gene2B),"Gene2B"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene2B),"Gene2B"])
rec[grep("\\(", rec$Gene1A),"Gene1A"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene1A),"Gene1A"])
rec[grep("\\(", rec$Gene2A),"Gene2A"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene2A),"Gene2A"])



rec4 <- rec %>% group_by(Sample, Gene1A,Gene1B) %>% summarise(Gene1.ct = n()) %>% filter(Gene1.ct >= 5) %>% as.data.frame() %>% cbind("multi_fused"="Gene1-multi_fused")
rec5 <- rec %>% group_by(Sample, Gene2A,Gene2B) %>% summarise(Gene2.ct = n()) %>% filter(Gene2.ct >= 5) %>% as.data.frame() %>% cbind("multi_fused"="Gene2-multi_fused")


rec4<-merge(total, rec4, by=c("Sample","Gene1A","Gene1B"))
rec5<-merge(total, rec5, by=c("Sample","Gene2A","Gene2B"))

colnames(rec4)
colnames(rec5)

total_multifused<-rbind(rec4[,c("Fused_Genes","Sample","Caller_type","Histology.Broad","note","multi_fused","JunctionReadCount","SpanningFragCount")],rec5[,c("Fused_Genes","Sample","Caller_type","Histology.Broad","note","multi_fused","JunctionReadCount","SpanningFragCount")]) 

total_multifused<-unique(total_multifused)




#filtered fusion through filtering strategy
saveRDS(total,"../../scratch/fusion_all_list_filt.rds")

#merge filtered fusion through literature and not multi_fused status
to.add$Caller_type<-paste(to.add$Caller,to.add$Fusion_Type,to.add$Confidence,sep="_")
final<-to.add[-which(to.add$Fused_Genes %in% total_multifused$Fused_Genes),]
final<-unique(final)


# save literature fused genes to run through expression filtering
saveRDS(final,"../../scratch/fusion_driver_list.rds")






