# This script outputs the oncoprint N counts table 

#### Set Up --------------------------------------------------------------------

# Load libraries
library(dplyr)
library(maftools)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to output directory for table produced
results_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "tables")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-m", "--maf_file_po"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains SNV information of primary_only samples",
  ),
  optparse::make_option(
    c("-c", "--cnv_file_po"),
    type = "character",
    default = NULL,
    help = "file path to file that contains CNV information of primary_only samples"
  ),
  optparse::make_option(
    c("-f", "--fusion_file_po"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information of primary_only samples"
  ),
  optparse::make_option(
    c("-g", "--maf_file_pp"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains SNV information of primary-plus samples",
  ),
  optparse::make_option(
    c("-d", "--cnv_file_pp"),
    type = "character",
    default = NULL,
    help = "file path to file that contains CNV information of primary-plus samples"
  ),
  optparse::make_option(
    c("-e", "--fusion_file_pp"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information of primary-plus samples"
  ),
  optparse::make_option(
    c("-s", "--metadata_file"),
    type = "character",
    default = NULL,
    help = "file path to the histologies file"
  ),
  optparse::make_option(
    c("-l", "--goi_list"),
    type = "character",
    default = NULL,
    help = "comma-separated list of genes of interest files that contain the
            genes to include on oncoprint"
  ),
  optparse::make_option(
    c("-b", "--broad_histology_list"),
    type = "character",
    default = NULL,
    help = "comma-separated list of `broad_histology` value to output numbers entering oncoprint"
  ),
  optparse::make_option(
    c("-o", "--short_name_list"),
    type = "character",
    default = NULL,
    help = "comma-separated list of file prefix for output tables"
  ),
  optparse::make_option(
    c("--include_introns"),
    action = "store_true",
    default = FALSE,
    help = "logical statement on whether to include intronic variants in oncoprint plot"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Define cnv_file, fusion_file, and genes object here as they still need to
# be defined for the `prepare_and_plot_oncoprint` custom function (the
# cnv_file specifically for the `read.maf` function within the custom function),
# even if they are NULL
broad_histology_list<-unlist(strsplit(opt$broad_histology_list,","))
short_name_list<-unlist(strsplit(opt$short_name_list,","))

# Source the custom functions script
source(
  file.path(
    root_dir,
    "analyses",
    "oncoprint-landscape",
    "util",
    "oncoplot-functions.R"
  )
)

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <- readr::read_tsv(opt$metadata_file, guess_max = 10000) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF file
maf_df_po <- data.table::fread(opt$maf_file_po,
                               stringsAsFactors = FALSE,
                               data.table = FALSE)
maf_df_pp <- data.table::fread(opt$maf_file_pp,
                               stringsAsFactors = FALSE,
                               data.table = FALSE)

if (!opt$include_introns) {
  maf_df_po <- maf_df_po %>%
    dplyr::filter(Variant_Classification != "Intron")
  
  maf_df_pp <- maf_df_pp %>%
    dplyr::filter(Variant_Classification != "Intron")
}

# Read in cnv file
cnv_df_po <- readr::read_tsv(opt$cnv_file_po) 
cnv_df_pp <- readr::read_tsv(opt$cnv_file_pp) 


# Read in fusion file and join
fusion_df_po <- readr::read_tsv(opt$fusion_file_po)
fusion_df_pp <- readr::read_tsv(opt$fusion_file_pp)


#### Set up oncoprint annotation objects --------------------------------------
oncoprint_n_table_po <- data.frame(matrix(ncol = 3, nrow = 0))
oncoprint_n_table_pp <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(oncoprint_n_table_po) <- c("broad_histology", "n_sample", "tumor_type")
colnames(oncoprint_n_table_pp) <- c("broad_histology", "n_sample", "tumor_type")

for(i in 1:length(broad_histology_list)){
  histology_of_interest <- broad_histology_list[i]
  
  if (!histology_of_interest == "Other CNS") {
    metadata_each <- metadata %>%
      dplyr::filter(broad_histology == histology_of_interest)
  } else {
    metadata_each <- metadata %>%
      dplyr::filter(
        broad_histology %in% c(
          "Ependymal tumor",
          "Tumors of sellar region",
          "Neuronal and mixed neuronal-glial tumor",
          "Tumor of cranial and paraspinal nerves",
          "Meningioma",
          "Mesenchymal non-meningothelial tumor",
          "Germ cell tumor",
          "Choroid plexus tumor",
          "Histiocytic tumor",
          "Tumor of pineal region",
          "Metastatic tumors",
          "Other astrocytic tumor",
          "Lymphoma",
          "Melanocytic tumor",
          "Other tumor"
        )
      )
    }
  # Now filter the remaining data files
  maf_po_each <- maf_df_po %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  cnv_po_each <- cnv_df_po %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  fusion_po_each <- fusion_df_po %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  # get bs ids in each dataframe
  maf_po_bsid <- maf_po_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  cnv_po_bsid <- cnv_po_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  fusion_po_bsid <- fusion_po_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  # take union
  all_po_bsid <- union(maf_po_bsid, cnv_po_bsid) %>% union(fusion_po_bsid) 
  # write out broad histology and number of samples to a table 
  oncoprint_n_table_po[i,1] <- histology_of_interest
  oncoprint_n_table_po[i,2] <- length(all_po_bsid)
  oncoprint_n_table_po[i,3] <- "primary_only"
  
  # Now do the same thing for primary plus samples
  maf_pp_each <- maf_df_pp %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  cnv_pp_each <- cnv_df_pp %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  fusion_pp_each <- fusion_df_pp %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata_each$Tumor_Sample_Barcode)
  
  # get bs ids in each dataframe
  maf_pp_bsid <- maf_pp_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  cnv_pp_bsid <- cnv_pp_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  fusion_pp_bsid <- fusion_pp_each %>% pull(Tumor_Sample_Barcode) %>% unique()
  # take union
  all_pp_bsid <- union(maf_pp_bsid, cnv_pp_bsid) %>% union(fusion_pp_bsid) 
  # write out broad histology and number of samples to a table 
  oncoprint_n_table_pp[i,1] <- histology_of_interest
  oncoprint_n_table_pp[i,2] <- length(all_pp_bsid)
  oncoprint_n_table_pp[i,3] <- "primary-plus"
  
}

# write out n number entering oncoprint
bind_rows(oncoprint_n_table_pp, oncoprint_n_table_po) %>% 
  readr::write_tsv(file.path(results_dir, "sample_n_in_oncoprint.tsv"))

