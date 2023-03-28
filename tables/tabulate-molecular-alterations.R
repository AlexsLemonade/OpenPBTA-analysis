# S. Spielman 2023 for CCDL
# 
# This script tabulates, for _all samples_ in OpenPBTA cohort, the presence/absence of 
#  alterations depicted in the manuscript's oncoprint plots.
# This script is heavily inspired by:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/oncoprint-landscape/01-map-to-sample_id.R
# 
# Output file structure:
#  sample_id    DNA BS     RNA BS     mapping_type       <------ alterations ------>
#  XXXX-YYYY    BS_###     BS_###    <1:1 or 1:many>            Y/N values
#
# Note that the oncoprint module is assumed to have been run before this script, 
#  such that `scratch/oncoprint-plots/` is populated.

# Define pipe
`%>%` <- dplyr::`%>%`

# Set up paths ------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
oncoprint_dir <- file.path(root_dir,
                           "analyses", 
                           "oncoprint-landscape")

pal_file           <- file.path(root_dir, "figures", "palettes", "broad_histology_cancer_group_palette.tsv")
metadata_file      <- file.path(data_dir, "pbta-histologies.tsv")
maf_file           <- file.path(data_dir, "pbta-snv-consensus-mutation.maf.tsv.gz")
hotspots_maf_file  <- file.path(data_dir, "pbta-snv-scavenged-hotspots.maf.tsv.gz")
cnv_autosomes_file <- file.path(data_dir, "consensus_seg_annotated_cn_autosomes.tsv.gz")
cnv_xy_file        <- file.path(data_dir, "consensus_seg_annotated_cn_x_and_y.tsv.gz")
fusion_file        <- file.path(data_dir, "pbta-fusion-putative-oncogenic.tsv")

  

#### Read in data ------------------------------------
histologies_df <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  dplyr::inner_join(
    readr::read_tsv(pal_file) %>%
      dplyr::select(cancer_group, cancer_group_display, broad_histology, broad_histology_display)
  )

# Only keep required maf columns 
keep_cols <-c("Hugo_Symbol", 
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2",
              "Variant_Classification",
              "Variant_Type",
              "Tumor_Sample_Barcode",
              "HGVSp_Short")

maf_df <- readr::read_tsv(maf_file) %>%
  dplyr::select(keep_cols)

hotspots_maf_df <- readr::read_tsv(hotspots_maf_file) %>%
  dplyr::select(keep_cols)
  
# merge hotspots maf to input maf
maf_df <- maf_df %>%
  dplyr::bind_rows(hotspots_maf_df) %>%
  dplyr::distinct()

# read and join cnv with histology information
cnv_autosomes_df <- readr::read_tsv(cnv_autosomes_file) %>%
  dplyr::left_join(
    dplyr::select(histologies_df, 
                  Kids_First_Biospecimen_ID,
                  germline_sex_estimate),
            by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")
  )

# combine autosomes with x/y cnv
cnv_df <- dplyr::bind_rows(
  cnv_autosomes_df, 
  readr::read_tsv(cnv_xy_file)
)

# read fusion data
fusion_df <- readr::read_tsv(fusion_file)



## Filter to relevant tumors --------------------------------------

# The final output should only tabuluate alterations for N=1074 tumors 
# We'll need to identify those samples, including their sample_id and biospecimen.

# In addition, we'll note if a sample_id is ambiguous. Ambiguous sample_id's will 
# have more than 2 rows associated with it in the histologies file when looking at 
# tumor samples -- that means we won't necessarily be able to determine **when an 
# WGS/WXS assay maps to an RNA-seq assay** for the purpose of tabulating molecular alterations

tumor_sample_ids_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue") %>%
  # keep these columns of interest moving forward:
  dplyr::select(sample_id, 
                Kids_First_Biospecimen_ID, 
                cancer_group_display, 
                broad_histology_display, 
                experimental_strategy) %>%
  dplyr::distinct() %>%
  # annotate if it is ambiguous or: any rows with >2 sample_ids
  # Associated with 2 sample_ids that have the same experimental strategy
  dplyr::add_count(sample_id) %>%
  dplyr::mutate(ambiguous = n > 2) %>% 
  # now, annotate if it has multiple rows of SAME experimental_strategy (needed for later wrangling)
  dplyr::add_count(sample_id, experimental_strategy) %>%
  dplyr::mutate(multiple_exp_strategy = n > 2) %>% 
  # remove temporary counting column
  dplyr::select(-n) %>%
  # add a modality column for later bookkeeping
  dplyr::mutate(modality = dplyr::case_when(
    experimental_strategy == "RNA-Seq" ~ "RNA",
    experimental_strategy == "WGS"     ~ "DNA",
    experimental_strategy == "WXS"     ~ "DNA",
    experimental_strategy == "Targeted Sequencing" ~ "DNA"
  )) %>%
  # remove no longer needed column
  dplyr::select(-experimental_strategy) 

# Pull out biospecimen id column for convenience 
bs_ids <- tumor_sample_ids_df %>%
  dplyr::pull(Kids_First_Biospecimen_ID) %>%
  unique() 

# Filter dfs to those ids:
histologies_df <- histologies_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% bs_ids)
maf_df <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% bs_ids)
cnv_df <- cnv_df %>%
  dplyr::filter(biospecimen_id %in% bs_ids)
fusion_df <- fusion_df %>%
  dplyr::filter(Sample %in% bs_ids)



#### MAF df preparation ------------------------------------------------

# join the sample_id information to the MAF file and then set as the tumor
# sample barcode
maf_df <- maf_df %>%
  dplyr::inner_join(
    tumor_sample_ids_df,
    by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")
  ) %>%
  # now let's remove this `Tumor_Sample_Barcode` column with biospecimen IDs in
  # preparation for our next step -- renaming `sample_id`
  dplyr::select(-Tumor_Sample_Barcode) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)



#### Fusion df preparation ---------------------------------------------------

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


### CNV df preparation -------------------------------------
cnv_df <- cnv_df %>%
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

### Combine processed data into one grand molecular alterations table ---------------------

# Now, we'll want to combine the following:
#   - maf_df
#   - reformat_fusion_df
#   - cnv_df
# As part of this, we'll keep the shared columns Tumor_Sample_Barcode, Hugo_Symbol, and Variant_Classification
# Then, we'll bring in tumor_sample_ids_df to match back to BS ids and diagnoses
# We'll have to process ambiguous and non-ambiguous separately and then merge back together

select_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")

maf_df <- maf_df %>%
  dplyr::select(select_cols)
reformat_fusion_df <- reformat_fusion_df %>%
  dplyr::select(select_cols)
cnv_df <- cnv_df %>%
  dplyr::select(select_cols)

# Bind all the rows:
alteration_df <- dplyr::bind_rows(maf_df, 
                                  reformat_fusion_df, 
                                  cnv_df) %>%
  # rename Tumor_Sample_Barcode to sample_id
  dplyr::rename(sample_id = Tumor_Sample_Barcode) %>%
  # testing
  dplyr::filter(sample_id %in% sample(unique(tumor_sample_ids_df$sample_id),10)) %>%
  dplyr::arrange(Hugo_Symbol)
#stop()

all_genes <- unique(alteration_df$Hugo_Symbol)
all_alts <- unique(alteration_df$Variant_Classification)
blank_alteration_df <- expand.grid(all_genes, all_alts) %>%
  purrr::set_names(c("Hugo_Symbol", "Variant_Classification")) %>%
  dplyr::mutate(present = "N") 


# Split into list of tibbles, one per sample_id
sample_alteration_df_list <- as.list(
  dplyr::group_split(alteration_df, sample_id)
) 

process_one_sample <- function(df, blank_alteration_df) {
  df %>%
    dplyr::mutate(present = "Y") %>%
    dplyr::full_join(blank_alteration_df) 
}







dplyr::mutate()
# This crashes R with the whole dataset. The operation doesn't use more than 8 GB RAM, 
# but somewhere it is getting upset.
# We'll need a diffferent (less elegant..) approach here to appease R.
# tidyr::complete(sample_id, Hugo_Symbol, Variant_Classification, 
#                 fill = list(present = "N")) %>%



# Determine whether each alteration is present in each sample_id
alteration_presence_df <- alteration_df %>%
  # If the row is here already, then the alteration is present
  dplyr::mutate(present = "Y") %>%
  dplyr::full_join(blank_alteration_df)



  # running this on the WHOLE thing crashes R. great stuff.
  # fill in rest of gene/alteration combinations, filling in present="N" (no) instead of `NA`
  
  # remove duplicate rows, since we care about presence/absense, not how many hits
  dplyr::distinct() %>%
  dplyr::ungroup() 

  
# enwiden: one column for each alteration type, with rows for samples/genes
  #tidyr::spread(Variant_Classification, present) %>%
  #


#############################################################################################
# Key result:
# `alteration_presence_df` holds a WIDE table of whether alterations are present/absent

#> colnames(alteration_presence_df)
# [1] "sample_id"              "Hugo_Symbol"           
# [3] "3'Flank"                "3'UTR"                 
# [5] "5'Flank"                "5'UTR"                 
# [7] "Amp"                    "Del"                   
# [9] "Frame_Shift_Del"        "Frame_Shift_Ins"       
#[11] "Fusion"                 "IGR"                   
#[13] "In_Frame_Del"           "In_Frame_Ins"          
#[15] "Intron"                 "Missense_Mutation"     
#[17] "Multi_Hit_Fusion"       "Nonsense_Mutation"     
#[19] "Nonstop_Mutation"       "RNA"                   
#[21] "Silent"                 "Splice_Region"         
#[23] "Splice_Site"            "Translation_Start_Site"
#############################################################################################


# Finally, we have to add/map in the biospecimen IDs. -----------------

# To begin, let's split out the fully identifiable samples
# For these, we can simply join in the ids
simple_ids <- tumor_sample_ids_df %>%
  dplyr::filter(!(ambiguous), 
                !(multiple_exp_strategy)) %>%
  dplyr::pull(sample_id) 

alteration_presence_df %>%
  dplyr::filter(sample_id %in% simple_ids) %>%
  dplyr::inner_join(tumor_sample_ids_df, 
                    by = "sample_id")


alteration_presence_df %>%
  dplyr::filter

alteration_presence_df %>%
  dplyr::ungroup() %>%
  dplyr::slice(116:121)

























# Combine with additional sample information: Fully identifiable samples
rna_ids <- tumor_sample_ids_df %>%
  dplyr::filter(!(ambiguous), 
                !(multiple_exp_strategy),
                modality == "RNA") %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID) 


dna_ids <- tumor_sample_ids_df %>%
  dplyr::filter(!(ambiguous), 
                !(multiple_exp_strategy),
                modality == "DNA") %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) 

  


