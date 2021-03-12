# K. S. Gaonkar 2021
# In this script, we will filter maf from a caller into 
# a resulting filtered file that has all calls that 
# - overlap with amino acid positions in a curated and 
# published cancer hotspot [database](https://www.cancerhotspots.org/files/hotspots_v2.xls)
# - overlap with non-coding region hotspots mutations ,


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list(
  make_option(c("-m", "--maffile"),type="character",
              help="Filtered calls from [mutect2|strelka2|vardict|lancet]"),
  make_option(c("-c", "--caller"), type="character",
              help="Caller type [mutect2|strelka2|vardict|lancet]"),
  make_option(c("-f","--cancer_hotspot_file"),
              help="MSKCC hotspot file", type="character"),
  make_option(c("-g","--genomic_site_hotspot_file"),
              help="genomic hotspot location file", type="character")
)

# Get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))
caller <- opt$caller
# Converting caller ID to call Uppercase for error handling
caller<-tolower(caller)

maf_coltypes = c(
  "Hugo_Symbol"="character",
  "Entrez_Gene_Id"="character",
  "Center" ="character",
  "NCBI_Build"="character",
  "Chromosome" ="character",
  "Start_Position" ="numeric",
  "End_Position" = "numeric",
  "Strand" = "character",
  "Variant_Classification" = "character",
  "Variant_Type" = "character",
  "Reference_Allele" = "character",
  "Tumor_Seq_Allele1" = "character",
  "Tumor_Seq_Allele2" = "character",
  "dbSNP_RS" = "character",
  "dbSNP_Val_Status" = "character",
  "Tumor_Sample_Barcode" = "character",
  "Matched_Norm_Sample_Barcode" = "character",
  "Match_Norm_Seq_Allele1" = "character",
  "Match_Norm_Seq_Allele2" = "character",
  "Tumor_Validation_Allele1" = "character",
  "Tumor_Validation_Allele2" = "character",
  "Match_Norm_Validation_Allele1" = "character",
  "Match_Norm_Validation_Allele2" = "character",
  "Verification_Status" = "character",
  "Validation_Status" = "character",
  "Mutation_Status" = "character",
  "Sequencing_Phase" = "character",
  "Sequence_Source" = "character",
  "Validation_Method" = "character",
  "Score" = "character",
  "BAM_File" = "character",
  "Sequencer" = "character",
  "Tumor_Sample_UUID" = "character",
  "Matched_Norm_Sample_UUID" = "character",
  "HGVSc" = "character",
  "HGVSp" = "character",
  "HGVSp_Short" = "character",
  "Transcript_ID" = "character",
  "Exon_Number" = "character",
  "t_depth" = "numeric",
  "t_ref_count" = "numeric",
  "t_alt_count" = "numeric",
  "n_depth" = "numeric",
  "n_ref_count" = "numeric",
  "n_alt_count" = "numeric",
  "all_effects" = "character",
  "Allele" = "character",
  "Gene" = "character",
  "Feature" = "character",
  "Feature_type" = "character",
  "Consequence" = "character",
  "cDNA_position" = "character",
  "CDS_position" = "character",
  "Protein_position" = "character",
  "Amino_acids" = "character",
  "Codons" = "character",
  "Existing_variation" = "character",
  "ALLELE_NUM" = "numeric",
  "DISTANCE" = "numeric",
  "STRAND_VEP" = "character",
  "SYMBOL" = "character",
  "SYMBOL_SOURCE" = "character",
  "HGNC_ID" = "character",
  "BIOTYPE" = "character",
  "CANONICAL" = "character",
  "CCDS" = "character",
  "ENSP" = "character",
  "SWISSPROT" = "character",
  "TREMBL" = "character",
  "UNIPARC" = "character",
  "RefSeq" = "character",
  "SIFT" = "character",
  "PolyPhen" = "character",
  "EXON" = "character",
  "INTRON" = "character",
  "DOMAINS" = "character",
  "AF" = "character",
  "AFR_AF" = "character",
  "AMR_AF" = "character",
  "ASN_AF" = "character",
  "EAS_AF" = "character",
  "EUR_AF" = "character",
  "SAS_AF" = "character",
  "AA_AF" = "character",
  "EA_AF" = "character",
  "CLIN_SIG" = "character",
  "SOMATIC" = "character",
  "PUBMED" = "character",
  "MOTIF_NAME" = "character",
  "MOTIF_POS" = "numeric",
  "HIGH_INF_POS" = "character",
  "MOTIF_SCORE_CHANGE"="double",
  "IMPACT" = "character",
  "PICK" = "character",
  "VARIANT_CLASS" = "character",
  "TSL" = "character",
  "HGVS_OFFSET" = "character",
  "PHENO" = "character",
  "MINIMISED" = "character",
  "ExAC_AF" = "character",
  "ExAC_AF_AFR" = "character",
  "ExAC_AF_AMR" = "character",
  "ExAC_AF_EAS" = "character",
  "ExAC_AF_FIN" = "character",
  "ExAC_AF_NFE" = "character",
  "ExAC_AF_OTH" = "character",
  "ExAC_AF_SAS" = "character",
  "GENE_PHENO" = "character",
  "FILTER" = "character",
  "flanking_bps" = "character",
  "vcf_id" = "character",
  "vcf_qual" = "double",
  "ExAC_AF_Adj" = "character",
  "ExAC_AC_AN_Adj" = "character",
  "ExAC_AC_AN" = "character",
  "ExAC_AC_AN_AFR" = "character",
  "ExAC_AC_AN_AMR" = "character",
  "ExAC_AC_AN_EAS" = "character",
  "ExAC_AC_AN_FIN" = "character",
  "ExAC_AC_AN_NFE" = "character",
  "ExAC_AC_AN_OTH" = "character",
  "ExAC_AC_AN_SAS" = "character",
  "ExAC_FILTER" = "character",
  "gnomAD_AF" = "character",
  "gnomAD_AFR_AF" = "character",
  "gnomAD_AMR_AF" = "character",
  "gnomAD_ASJ_AF" = "character",
  "gnomAD_EAS_AF" = "character",
  "gnomAD_FIN_AF" = "character",
  "gnomAD_NFE_AF" = "character",
  "gnomAD_OTH_AF" = "character",
  "gnomAD_SAS_AF" = "character",
  "vcf_pos" = "numeric"
  )

maffile <- data.table::fread(opt$maffile,
                             stringsAsFactors = FALSE,
                             col.names = names(maf_coltypes),
                             colClasses = maf_coltypes)

mskcc_cancer_hotspot_file <- opt$cancer_hotspot_file
genomic_region_file <- opt$genomic_site_hotspot_file


# Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir,
                         "analyses",
                         "hotspots-detection",
                         "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# MSKCC cancer hotspot file is divided into snv and indels
# in 2 tabs in the excel file provided in their [database](https://www.cancerhotspots.org/#/download) 
# We will combine the 2 sheets to create a `hotspot_database_amino_acid` dataframe 

# MSKCC cancer database
hotspot_database_2017_snv <- readxl::read_xls(mskcc_cancer_hotspot_file,sheet = 1) 
hotspot_database_2017_indel <- readxl::read_xls(mskcc_cancer_hotspot_file,sheet = 2)
hotspot_database_amino_acid <- bind_rows( hotspot_database_2017_snv, hotspot_database_2017_indel) %>%
  dplyr::select("Amino_Acid_Position","Hugo_Symbol") %>%
  dplyr::mutate(hotspot_database="MSKCC") %>%
  unique() %>%
  as.data.frame()

# TERT promoter region
hotspot_database_genomic <- read_tsv(genomic_region_file)

# Gather Hugo_Symnol from both cancer hotspot and genomic regions to filter maf files
# This way we can directly filter large mafs into smaller manageable sizes to 
# check for overlaps
gene_table <- data.frame("Hugo_Symbol"=c(hotspot_database_amino_acid$Hugo_Symbol,
                                         hotspot_database_genomic$Hugo_Symbol))

print(paste("Subsetting" ,caller, "maf for hotspots"))

source ("utils/prepMaf.R")
# prepMaf will filter the input given maf file, given a
# - dataframe of hotspot_amino_acid_position_df to 
# check for Amino_Acid_Position overlaps per Hugo_Symbol
# - dataframe of hotspot_genomic_site_df to check for
# region overlap in the given Start_Position and End_Position
# - impact_values to filter IMPACT column in maf file
# - gene_table to filter Hugo_Symbol in maf file

maf_subset<-prepMaf(maffile, gene_table = gene_table,
                        impact_values = "MODERATE|HIGH|MODIFIER",
                        hotspot_amino_acid_position_df = hotspot_database_amino_acid, 
                        hotspot_genomic_site_df = hotspot_database_genomic) %>%
  mutate(caller=caller)	

# save 
saveRDS(maf_subset, file.path(results_dir,paste0(caller,"_hotspots.RDS")))
