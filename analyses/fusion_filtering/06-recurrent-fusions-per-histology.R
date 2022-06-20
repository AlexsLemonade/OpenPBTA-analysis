# K. S. Gaonkar 2019
# Identify recurrent fusion and genes per broad histology
# 
# Sample selection criteria : removed cell-lines to only keep tumor samples 

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))

option_list <- list(
  make_option(c("-n", "--seed"),type = "integer",
             help="seed integer",default = 2020),
  make_option(c("-S", "--standardFusionCalls"),type="character",
              help="Standardized fusion calls (.tsv) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-o","--outputfolder"),type="character",
              help="output folder for results "),
  make_option(c("-i","--independentSpecimensFile"),type="character",
              help="independent specimens WGS/WXS to match to rnaseq")
)

opt <- parse_args(OptionParser(option_list=option_list))
standardFusionCalls<-opt$standardFusionCalls
clinicalFile<-opt$clinicalFile
outputfolder<-opt$outputfolder
independentSpecimensFile<-opt$independentSpecimensFile
seed <- opt$seed

set.seed(seed)

# input data
standardFusionCalls <- read_tsv(standardFusionCalls) %>% dplyr::arrange(Sample) %>% as.data.frame()
clinical<-read_tsv(clinicalFile,col_types = readr::cols(molecular_subtype = readr::col_character())) %>% arrange(Kids_First_Biospecimen_ID,sample_id,tumor_descriptor,experimental_strategy,composition)
# gather RNA-seq from WGS/WXS samples in independent-specimens.wgswxs.primary-plus.tsv
independentSpecimens<-read_tsv(independentSpecimensFile) %>% arrange(Kids_First_Biospecimen_ID) %>% as.data.frame()

sampleIDMatchedIndependent<-clinical %>% dplyr::filter(Kids_First_Biospecimen_ID %in% independentSpecimens$Kids_First_Biospecimen_ID) %>% dplyr::select(sample_id) %>% as.vector() 

# PNOC sampleIDs have .WXS which would need to .RNA-Seq ; Panel not in independent sample list ; so using patient ID to match
clinical_pnoc<-clinical %>% 
  dplyr::filter(cohort=="PNOC",experimental_strategy=="RNA-Seq",tumor_descriptor=="Initial CNS Tumor") %>% 
  dplyr::select(Kids_First_Participant_ID,Kids_First_Biospecimen_ID)


clinical_rna<-clinical %>% 
  # select WGS/WXS matched sampleIDs + experimental_strategy=="RNA-Seq" +composition == "Solid Tissue"
  dplyr::filter(sample_id %in% sampleIDMatchedIndependent$sample_id,
         experimental_strategy == "RNA-Seq",
         composition == "Solid Tissue") %>%
  dplyr::group_by(Kids_First_Participant_ID) %>%
  # sample 1 of BS_IDs for multiple sampleID found in matching with WGS/WXS
  dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1)) 

# match RNASeq only files
clinical_wgs<-clinical %>% 
  dplyr::filter(experimental_strategy == "WGS" | experimental_strategy == "WXS",sample_type=="Tumor")

# remove samples which have WGS/WXS because those would have been captured from the independent-wgswxs-sample
clinical_rna_v2<-clinical %>% 
  dplyr::filter(experimental_strategy == "RNA-Seq",!Kids_First_Participant_ID %in% clinical_wgs$Kids_First_Participant_ID)

clinical_rna_intial<-clinical_rna_v2 %>% 
  dplyr::filter(composition=="Solid Tissue",tumor_descriptor=="Initial CNS Tumor") %>% 
  dplyr::select("Kids_First_Participant_ID","Kids_First_Biospecimen_ID")

clinical_rna_non_initial<- clinical_rna_v2 %>% 
  dplyr::filter(!Kids_First_Participant_ID %in% clinical_rna_intial$Kids_First_Participant_ID ) %>% 
  dplyr::filter(experimental_strategy=="RNA-Seq" ,composition=="Solid Tissue") %>%
  as.data.frame() %>% 
  dplyr::group_by(Kids_First_Participant_ID) %>%
  # random sample 1 of BS_IDs for multiple sampleID found in matching with WGS/WXS
  dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1)) 

# bind WGS/WXS matched RNA-Seq + RNAseq only samples initial tumor samples + RNAseq only recurrent/progressive samples + RNASeq matched to PNOC samples
clinical_rna<-rbind(clinical_rna,clinical_rna_intial,clinical_rna_non_initial,clinical_pnoc) %>% unique()

clinical_rna<-clinical_rna %>% left_join(clinical,by=c("Kids_First_Participant_ID","Kids_First_Biospecimen_ID"))


# Putative Driver Fusions annotated with broad_histology
standardFusionCalls<-standardFusionCalls %>% 
  dplyr::filter( Sample%in% clinical_rna$Kids_First_Biospecimen_ID) %>%
  left_join(clinical,by=c("Sample"="Kids_First_Biospecimen_ID","Kids_First_Participant_ID")) %>% 
  dplyr::filter(!is.na(broad_histology)) %>% as.data.frame()

# keep only inframe fusions 
standardFusionCalls <- standardFusionCalls %>%
  # filter to keep inframe fusions 
  dplyr::filter(Fusion_Type %in% c("in-frame"))

# running this to remove GeneA==GeneB which might be non-canonical transcripts/ false positives from arriba documentation ; still unknown function/relevance
standardFusionCalls<-standardFusionCalls %>%
  dplyr::mutate(removal = dplyr::case_when(
    Gene1A == Gene2A | Gene1A == Gene2B | Gene2A == Gene1B | Gene1A == Gene1B | Gene2A == Gene2B ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  # remove the ones you want to remove by negating removal
  dplyr::filter(!removal) %>%
  # get rid of removal column
  dplyr::select(-removal)

# gather recurrent fusion per patient per broad_histology
rec_fusions <- standardFusionCalls %>% 
  dplyr::select("FusionName","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  dplyr::group_by(FusionName, broad_histology) %>% 
  dplyr::summarize(count = n())

#find rec fusions per PATIENT per broad_histology
rec_fusions<-rec_fusions[rec_fusions$count>3,]
rec_fusions<-rec_fusions[order(rec_fusions$count,decreasing = TRUE),]
write.table(rec_fusions,file.path(outputfolder,"pbta-fusion-recurrent-fusion-byhistology.tsv"),quote = FALSE,row.names = FALSE,sep="\t")


# binary matrix for recurrent fusions found in SAMPLE per broad_histology
rec_fusions_mat<-rec_fusions %>% 
  # to add sample with rec fusion
  left_join(standardFusionCalls, by=c("FusionName","broad_histology")) %>%
  select("FusionName","broad_histology","Sample") %>% 
  # to add all rna Sample for binary matrix
  full_join(clinical_rna,by=c("Sample"="Kids_First_Biospecimen_ID","broad_histology"="broad_histology"))

# adding a full join to recurrent fusion adds NAs to FusionName column that don't have recurrent fusions. Changing NA to No_rec_fusions so in matrix format it columns specifies that instead of NA
rec_fusions_mat[is.na(rec_fusions_mat$FusionName),"FusionName"]<-"No_rec_fusion"

# binary matrix
rec_fusions_mat<-dcast(rec_fusions_mat,Sample~FusionName,value.var = "Sample",fun.aggregate = function(x){as.integer(length(x) > 0)},drop = FALSE) 

write.table(rec_fusions_mat,file.path(outputfolder,"pbta-fusion-recurrent-fusion-bysample.tsv"),quote = FALSE,row.names = FALSE,sep="\t")


# gene1A recurrent
rec_gene1A <- standardFusionCalls %>% 
  dplyr::select("Gene1A","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  dplyr::rename("Gene"="Gene1A")

# gene1B recurrent
rec_gene1B <- standardFusionCalls %>% 
  dplyr::select("Gene1B","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  dplyr::rename("Gene"="Gene1B")

# merge and then summarize
rec_gene<-rbind(rec_gene1A,rec_gene1B) %>% 
  unique() %>%
  dplyr::group_by(Gene,broad_histology) %>% 
  dplyr::summarize(count = n()) 

#find rec fused genes per PATIENT per broad_histology
rec_gene<-rec_gene[rec_gene$count>3,]
rec_gene<-rec_gene[order(rec_gene$count,decreasing = TRUE),]
write.table(rec_gene,file.path(outputfolder,"pbta-fusion-recurrently-fused-genes-byhistology.tsv"),quote = FALSE,row.names = FALSE,sep="\t")

# binary matrix for recurrently fused genes found in SAMPLE per broad_histology
rec_geneA_mat<-rec_gene %>% 
  # to add sample with rec fusion gene A (5' gene)
  left_join(standardFusionCalls, by=c("Gene"="Gene1A","broad_histology")) %>%
  select("Gene","broad_histology","Sample") %>%
  filter(!is.na(Sample))

rec_geneB_mat<-rec_gene %>% 
  # to add sample with rec fusion gene B (3' gene)
  left_join(standardFusionCalls, by=c("Gene"="Gene1B","broad_histology")) %>%
  select("Gene","broad_histology","Sample") %>%
  filter(!is.na(Sample))

rec_gene_mat<-rbind(rec_geneA_mat,rec_geneB_mat) %>% unique() %>%
  # to add all rna Sample for binary matrix
  full_join(clinical_rna,by=c("Sample"="Kids_First_Biospecimen_ID","broad_histology"="broad_histology")) %>%
  unique()
rec_gene_mat[is.na(rec_gene_mat$Gene),"Gene"]<-"No_rec_fused_gene"

# binary matrix
rec_gene_mat<-dcast(rec_gene_mat,Sample~Gene,value.var = "Sample",fun.aggregate = function(x){as.integer(length(x) > 0)},drop = FALSE)

write.table(rec_gene_mat,file.path(outputfolder,"pbta-fusion-recurrently-fused-genes-bysample.tsv"),quote = FALSE,row.names = FALSE,sep="\t")


