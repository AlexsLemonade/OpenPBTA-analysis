# K S Gaonkar 2020
# Subsets the consensus mutation file to biospecimens where short_histology == LGAT

# Look for git root folder
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# get subset folder
subset_dir<-file.path(root_dir,"analyses","molecular-subtyping-LGAT","lgat-subset")

# clinical file
clinical<- readr::read_tsv(file.path(root_dir,"data","pbta-histologies.tsv"))
# consensus mutation data
consensusMutation <- readr::read_tsv(file.path(root_dir,"data","pbta-snv-consensus-mutation.maf.tsv.gz"))


# Filter for LGAT samples
lgat_wgs_subset <- clinical %>%
  dplyr::filter(short_histology == "LGAT",
                sample_type == "Tumor",
                composition == "Solid Tissue",
                experimental_strategy == "WGS")

# Filter consensus mutation files for LGAT subset
consensusMutationSubset <- consensusMutation %>%
  # find lgat samples
  dplyr::filter(Tumor_Sample_Barcode %in% lgat_subset$Kids_First_Biospecimen_ID) %>%
  # look for BRAF V600E mutations and make a column for BRAF_V600E
  dplyr::filter(Hugo_Symbol == "BRAF" & HGVSp_Short == "p.V600E") %>%
  # select tumor sample barcode
  dplyr::select(Tumor_Sample_Barcode,HGVSp_Short) %>% 
  unique() %>%
  # join other WGS LGAT samples
  full_join(lgat_wgs_subset,by=c("Tumor_Sample_Barcode"="Kids_First_Biospecimen_ID")) %>%
  # get BRAF_V600E status
  dplyr::mutate(BRAF_V600E=dplyr::case_when(HGVSp_Short=="p.V600E"~"Yes",is.na(HGVSp_Short) ~"No"))

# remove consensusMutation
rm(consensusMutation)

# save to subset folder
readr::write_tsv(consensusMutationSubset,file.path(subset_dir, "LGAT_snv_subset.tsv"))
