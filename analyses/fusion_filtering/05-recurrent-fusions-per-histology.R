# K. S. Gaonkar 2019
# Identify recurrent fusion and genes per broad histology
# 
# Sample selection criteria : removed cell-lines to only keep tumor samples 

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))

option_list <- list(
  make_option(c("-S", "--standardFusionCalls"),type="character",
              help="Standardized fusion calls (.tsv) "),
  make_option(c("-c","--clinicalFile"),type="character",
              help="clinical file for all samples (.tsv)"),
  make_option(c("-o","--outputfolder"),type="character",
              help="output folder for results ")
)

opt <- parse_args(OptionParser(option_list=option_list))
standardFusionCalls<-opt$standardFusionCalls
clinicalFile<-opt$clinicalFile
outputfolder<-opt$outputfolder


# input data
standardFusionCalls <- read_tsv(standardFusionCalls) %>% as.data.frame()
clinical<-read_tsv(clinicalFile)

# to get initial tumor samples 
clinical_rna<-clinical %>% dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  # samples 
  dplyr::select("Kids_First_Biospecimen_ID","Kids_First_Participant_ID","source_text_tumor_descriptor","primary_site","composition","broad_histology") %>% unique() %>% 
  group_by(Kids_First_Participant_ID) %>%
  # filter cell-lines 
  dplyr::filter(composition=="Solid Tissue",source_text_tumor_descriptor=="Initial CNS Tumor") %>%
  dplyr::mutate(sample.count.per.pt_ID = n()) %>% as.data.frame()


#get other patients where initial tumors were not found
non_initial_tumor_samples<- clinical[-which(clinical$Kids_First_Participant_ID %in% clinical_rna$Kids_First_Participant_ID ),] %>% dplyr::filter(experimental_strategy=="RNA-Seq" ,composition=="Solid Tissue") %>% as.data.frame()


clinical_rna<-rbind(clinical_rna ,clinical %>% 
                      dplyr::filter(Kids_First_Biospecimen_ID %in% non_initial_tumor_samples$Kids_First_Biospecimen_ID) %>% 
                      # samples 
                      dplyr::select("Kids_First_Biospecimen_ID","Kids_First_Participant_ID","source_text_tumor_descriptor","primary_site","composition","broad_histology") %>% 
                      unique() %>% 
                      group_by(Kids_First_Participant_ID) %>%
                      dplyr::mutate(sample.count.per.pt_ID = n()) %>% 
                      as.data.frame() )
                    

# get samples with multiple tumors to review
clinical_rna<-clinical_rna[order(clinical_rna$Kids_First_Participant_ID,decreasing = TRUE),]
write.table(clinical_rna[clinical_rna$sample.count.per.pt_ID>1,],file.path(outputfolder,"multiple_tumor_rnaseq.tsv"),sep="\t",quote = FALSE,row.names = FALSE)

# Putative Driver Fusions annotated with broad_histology
standardFusionCalls<-standardFusionCalls %>% left_join(clinical_rna,by=c("Sample"="Kids_First_Biospecimen_ID","Kids_First_Participant_ID")) %>% dplyr::filter(!is.na(broad_histology)) %>% as.data.frame()

#remove fusions found in benign tumors as internal false positive control
# KCNH1--AL590132.1 and AL590132.1--KCNH1 come up as recurrent in multiple histologies but using the following filter helps remove such fusions.
standardFusionCalls<-standardFusionCalls[-which(standardFusionCalls$FusionName %in%
                                    unique(standardFusionCalls[which(standardFusionCalls$broad_histology=="Benign tumor"),"FusionName"])),] 


# only inframe 
standardFusionCalls<-standardFusionCalls %>% dplyr::filter(Fusion_Type %in% c("in-frame"))

# running this to remove GeneA==GeneB which might be non-canonical transcripts/ false positives from arriba documentation 
standardFusionCalls<-standardFusionCalls[-which(standardFusionCalls$Gene1A==standardFusionCalls$Gene2A|standardFusionCalls$Gene1A==standardFusionCalls$Gene2B|standardFusionCalls$Gene2A==standardFusionCalls$Gene1B|standardFusionCalls$Gene1A==standardFusionCalls$Gene1B|standardFusionCalls$Gene2B==standardFusionCalls$Gene2A),]


# gather recurrent fusion per patient per broad_histology
rec_fusions<-standardFusionCalls %>%
  as.data.frame() %>%
  dplyr::select("FusionName","broad_histology","Kids_First_Participant_ID") %>%
  unique() %>%
  group_by(FusionName,broad_histology) %>%
  dplyr::select(-Kids_First_Participant_ID) %>%
  mutate(count=n()) %>%
  unique() %>%
  as.data.frame()

#find rec fusions per PATIENT per broad_histology
rec_fusions<-rec_fusions[rec_fusions$count>3,]
rec_fusions_all<-rec_fusions[order(rec_fusions$count,decreasing = TRUE),]
write.table(rec_fusions_all,file.path(outputfolder,"rec_fusions_participant_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")

# binary matrix for recurrent fusions found in SAMPLE per broad_histology
rec_fusions_mat<-rec_fusions %>% 
  # to add sample with rec fusion
  left_join(standardFusionCalls, by=c("FusionName","broad_histology")) %>% select("FusionName","broad_histology","Sample") %>% 
  # to add all rna Sample for binary matrix
  full_join(clinical_rna,by=c("Sample"="Kids_First_Biospecimen_ID","broad_histology"="broad_histology"))

rec_fusions_mat[is.na(rec_fusions_mat$FusionName),"FusionName"]<-"No_rec_fusion"

# binary matrix
rec_fusions_mat<-dcast(rec_fusions_mat,Sample~FusionName,value.var = "Sample",fun.aggregate = function(x){as.integer(length(x) > 0)},drop = FALSE) 

write.table(rec_fusions_mat,file.path(outputfolder,"rec_fusions_matrix_sample_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")


# gene1A recurrent
rec_gene1A<-standardFusionCalls %>%
  dplyr::select("Gene1A","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  group_by(Gene1A,broad_histology) %>%
  dplyr::select(-Kids_First_Participant_ID) %>%
  mutate(count=n())%>% unique() %>% as.data.frame()

colnames(rec_gene1A)<-c("Gene","broad_histology","count")

# gene1B recurrent
rec_gene1B<-standardFusionCalls %>%
  dplyr::select("Gene1B","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  group_by(Gene1B,broad_histology) %>%
  dplyr::select(-Kids_First_Participant_ID) %>%
  mutate(count=n())%>% unique() %>% as.data.frame()

colnames(rec_gene1B)<-c("Gene","broad_histology","count")
rec_gene<-rbind(rec_gene1A,rec_gene1B)

#find rec fused genes per PATIENT per broad_histology
rec_gene<-rec_gene[rec_gene$count>3,]
rec_gene_all<-rec_gene[order(rec_gene$count,decreasing = TRUE),]
write.table(rec_gene,file.path(outputfolder,"rec_genes_participant_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")

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

write.table(rec_gene_mat,file.path(outputfolder,"rec_genes_matrix_sample_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")






