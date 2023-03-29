## S. Spielman 2023 from CCDL
## This file holds functions to prepare maf, cnv, and fusion data for Zenodo upload
## These functions are NOT TO BE USED outside of this context.
# This code is heavily/entirely inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R




prepare_maf <- function(maf_df) {
   maf_df %>%
    dplyr::inner_join(
      tumor_sample_ids_df,
      by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")
    ) %>%
    # now let's remove this `Tumor_Sample_Barcode` column with biospecimen IDs in
    # preparation for our next step -- renaming `sample_id`
    dplyr::select(-Tumor_Sample_Barcode) %>%
    dplyr::rename(Tumor_Sample_Barcode = sample_id)
}



prepare_fusion <- function(fusion_df) {
  
  # We'll handle fusions where reciprocal fusions exist (e.g., 
  # reciprocal_exists == TRUE) separately from other fusions
  fusion_reciprocal_df <- fusion_df %>%
    # Reciprocal fusions only
    dplyr::filter(reciprocal_exists) %>%
    # BSID + Gene1--Gene2
    dplyr::select(Sample, FusionName) %>%
    # Because we're only looking at presence or absence here, we can filter to 
    # distinct identifier-fusion pairs
    dplyr::distinct() %>%
    dplyr::group_by(Sample, FusionName) %>%
    # Put genes in the fusions in alphabetical order to collapse
    dplyr::mutate(SortedFusionName = stringr::str_c(
      sort( stringr::str_split(FusionName, "--", simplify = TRUE) ), 
      collapse = "--")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Sample,
                  FusionName = SortedFusionName) %>%
    # When fusions are in alphabetical order, this ensures each reciprocal
    # fusion is only counted once
    dplyr::distinct()
  
  # No need to put fusions where no reciprocal exists in alphabetical order,
  # so we handle these separately
  fusion_no_reciprocal_df <- fusion_df %>%
    dplyr::filter(!reciprocal_exists) %>%
    dplyr::select(Sample, FusionName) %>%
    # But we can remove duplicates to avoid them being counted as multihit
    dplyr::distinct()
  
  # Bind reciprocal and no reciprocal together
  fusion_filtered_df <- dplyr::bind_rows(fusion_reciprocal_df, fusion_no_reciprocal_df)
  rm(fusion_reciprocal_df, fusion_no_reciprocal_df)
  
  # Separate fusion gene partners
  fus_sep <- fusion_filtered_df %>%
    # Separate the 5' and 3' genes
    tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--") %>%
    # Use row numbers to mark unique fusions - this will help us when
    # we melt and remove selfie fusions below
    tibble::rowid_to_column("Fusion_ID") %>%
    dplyr::select(Fusion_ID, Sample, Gene1, Gene2) %>%
    reshape2::melt(id.vars = c("Fusion_ID", "Sample"),
                   variable.name = "Partner",
                   value.name = "Hugo_Symbol") %>%
    dplyr::arrange(Fusion_ID)
  
  multihit_fusions <- fus_sep %>%
    # For looking at multi-hit fusions, we want to remove selfie fusions
    # If we drop the 5', 3' information in Partner, we can use distinct
    # to remove the selfie fusions because Fusion_ID marks what fusion
    # it came from
    dplyr::select(-Partner) %>%
    dplyr::distinct() %>%
    # Now we want to count how many times a gene is fused within a sample
    # Anything with more than one count is considered multi-hit
    dplyr::count(Sample, Hugo_Symbol) %>%
    dplyr::filter(n > 1) %>%
    dplyr::select(-n) %>%
    dplyr::mutate(Variant_Classification = "Multi_Hit_Fusion")
  
  # Filter out multi-hit fusions from the other fusions that we will
  # label based on whether they are the 5' or 3' gene
  single_fusion <- fus_sep %>%
    dplyr::anti_join(multihit_fusions,
                     by = c("Sample", "Hugo_Symbol")) %>%
    dplyr::select(Sample, Hugo_Symbol) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Variant_Classification = "Fusion")
  
  # Combine multi-hit and single fusions into one data frame
  # And add the other identifiers!
  reformat_fusion_df <- multihit_fusions %>%
    dplyr::bind_rows(single_fusion) %>%
    dplyr::mutate(Variant_Type = "OTHER") %>%
    dplyr::inner_join(tumor_sample_ids_df,
                      by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
    dplyr::rename(Tumor_Sample_Barcode = sample_id,
                  Kids_First_Biospecimen_ID = Sample)

  # return reformatted
  reformat_fusion_df
}



prepare_cnv <- function(cnv_df) {
  
  cnv_df %>%
    dplyr::inner_join(
      tumor_sample_ids_df,
      by = c("biospecimen_id"="Kids_First_Biospecimen_ID")) %>%
    dplyr::mutate(Tumor_Sample_Barcode =  sample_id) %>%
    dplyr::rename(Variant_Classification = status,
                  Hugo_Symbol = gene_symbol) %>%
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
    # mutate loss and amplification to Del and Amp to fit Maftools format
    dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "deep deletion" ~ "Del",
                                                            Variant_Classification == "amplification" ~ "Amp",
                                                            TRUE ~ as.character(Variant_Classification))) %>%
    # only keep Del and Amp calls
    dplyr::filter(Variant_Classification %in% c("Del", "Amp"))
  
}
