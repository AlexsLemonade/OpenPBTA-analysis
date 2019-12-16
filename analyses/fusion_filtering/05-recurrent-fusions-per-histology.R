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

# filter to the tumor RNA-seq samples first 
clinical_rna <- clinical %>%
  dplyr::filter(experimental_strategy == "RNA-Seq",
                composition == "Solid Tissue") %>%
  # only the columns we will use downstream
  dplyr::select("Kids_First_Biospecimen_ID",
                "Kids_First_Participant_ID",
                "tumor_descriptor",
                "primary_site",
                "composition",
                "broad_histology")

# filter the tumor RNA-seq samples to only initial tumors
initial_tumor <- clinical_rna %>%
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor")

# use initial tumor participant IDs to filter to only 
# participants without initial tumors
non_initial_tumor <- clinical_rna %>%
  dplyr::filter(!(Kids_First_Participant_ID %in% initial_tumor$Kids_First_Participant_ID))

# join rows together for initial and non-initial
clinical_rna <- initial_tumor %>%
  dplyr::bind_rows(non_initial_tumor) %>%
  # count how many times a participant ID appears
  dplyr::group_by(Kids_First_Participant_ID) %>%
  dplyr::mutate(sample.count.per.pt_ID = n()) %>% 
  as.data.frame() %>%
  # order rows by participant ID
  dplyr::arrange(dplyr::desc(Kids_First_Participant_ID))
                    

# get samples with multiple tumors to review
clinical_rna<-clinical_rna[order(clinical_rna$Kids_First_Participant_ID,decreasing = TRUE),]
write.table(clinical_rna[clinical_rna$sample.count.per.pt_ID>1,],file.path(outputfolder,"multiple_tumor_rnaseq.tsv"),sep="\t",quote = FALSE,row.names = FALSE)

# Putative Driver Fusions annotated with broad_histology
standardFusionCalls<-standardFusionCalls %>% left_join(clinical_rna,by=c("Sample"="Kids_First_Biospecimen_ID","Kids_First_Participant_ID")) %>% dplyr::filter(!is.na(broad_histology)) %>% as.data.frame()

#remove fusions found in benign tumors as internal false positive control
# KCNH1--AL590132.1 and AL590132.1--KCNH1 come up as recurrent in multiple histologies but using the following filter helps remove such fusions.

# get the names of fusions that are present in benign tumors
fusions_in_benign <- standardFusionCalls %>%
  dplyr::filter(broad_histology == "Benign tumor") %>%
  unique() %>%
  dplyr::pull(FusionName)

# remove fusions found in benign tumors as internal false positive control
standardFusionCalls <- standardFusionCalls %>%
  dplyr::filter(!(FusionName %in% fusions_in_benign)) %>% 
  # keep only inframe fusions 
  dplyr::filter(Fusion_Type %in% c("in-frame"))

# running this to remove GeneA==GeneB which might be non-canonical transcripts/ false positives from arriba documentation 
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
write.table(rec_fusions,file.path(outputfolder,"rec_fusions_participant_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")

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

write.table(rec_fusions_mat,file.path(outputfolder,"rec_fusions_matrix_sample_histology_level.tsv"),quote = FALSE,row.names = FALSE,sep="\t")


# gene1A recurrent
rec_gene1A <- standardFusionCalls %>% 
  dplyr::select("Gene1A","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  dplyr::group_by(Gene1A,broad_histology) %>% 
  dplyr::summarize(count = n()) %>%
  dplyr::rename("Gene"="Gene1A")

# gene1B recurrent
rec_gene1B <- standardFusionCalls %>% 
  dplyr::select("Gene1B","broad_histology","Kids_First_Participant_ID")%>%
  unique() %>%
  dplyr::group_by(Gene1B,broad_histology) %>% 
  dplyr::summarize(count = n()) %>%
  dplyr::rename("Gene"="Gene1B")

rec_gene<-rbind(rec_gene1A,rec_gene1B)

#find rec fused genes per PATIENT per broad_histology
rec_gene<-rec_gene[rec_gene$count>3,]
rec_gene<-rec_gene[order(rec_gene$count,decreasing = TRUE),]
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





