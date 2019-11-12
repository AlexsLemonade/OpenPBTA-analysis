
# takes:
# cnv file
# fusion file
# snv file
# metadata file
# optionally, the independent specimens list

library(dplyr)

#### Command line options ------------------------------------------------------

option_list <- list(
  optparse::make_option(
    c("--maf_file"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains snv information",
  ),
  optparse::make_option(
    c("--cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to SEG file that contains cnv information"
  ),
  optparse::make_option(
    c("--fusion_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information"
  ),
  optparse::make_option(
    c("--metadata_file"),
    type = "character",
    default = NULL,
    help = "file path to clinical file"
  ),
  optparse::make_option(
    c("--filename_lead"),
    type = "character",
    default = NULL,
    help = ""
  ),
  optparse::make_option(
    c("--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory for files that are prepared in this script"
  ),
  optparse::make_option(
    c("--independent_specimens"),
    type = "character",
    default = NULL,
    help = "optional file path to independent specimens list"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### Directory set up ----------------------------------------------------------

output_dir <- opt$output_directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#### Read in data --------------------------------------------------------------

histologies_df <- readr::read_tsv(opt$metadata_file)
maf_df <- readr::read_tsv(opt$maf_file)

biospecimens_to_remove <- unique(c(problem_biospecimens,
                                   not_tumor_biospecimens))

maf_df <- maf_df %>%
  dplyr::filter(!(Tumor_Sample_Barcode %in% biospecimens_to_remove))

cnv_df <- readr::read_tsv(opt$cnv_file)
fusion_df <- readr::read_tsv(opt$fusion_file)

#### Filter to independent specimens (optional) --------------------------------

if (!is.null(opt$independent_specimens)) {

  independent_df <- readr::read_tsv(opt$independent_specimens)
  ind_biospecimen <- independent_df %>%
    pull(Kids_First_Biospecimen_ID)

  # filter the genome data, e.g., SNV and CNV data, to only include biospecimen
  # identifiers in the independent file
  maf_df <- maf_df %>%
    filter(Tumor_Sample_Barcode %in% ind_biospecimen)
  cnv_df <- cnv_df %>%
    filter(biospecimen_id %in% ind_biospecimen)

}

#### MAF file preparation ------------------------------------------------------

# we need to add the participant ID and sample_id, we will paste these together
# to use as the Tumor_Sample_Barcode
maf_df <- maf_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  rename(Tumor_Sample_Barcode = sample_id)

# Write MAF to file
maf_output <- file.path(output_dir, paste0(opt$filename_lead, "_maf.tsv"))
readr::write_tsv(maf_df, maf_output)

#### Fusion file preparation ---------------------------------------------------

# Separate fusion gene partners and add variant classification and center
fus_sep <- fusion_df %>%
  # Separate the 5' and 3' genes
  tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--") %>%
  select(Sample, Gene1, Gene2)

reformat_fusion <- fus_sep %>%
  # Here we want to tally how many times the 5' gene shows up as a fusion hit
  # in a sample
  group_by(Sample, Gene1) %>%
  tally() %>%
  # If the sample-5' gene pair shows up more than once, call it a multi hit
  # fusion
  mutate(Variant_Classification =
           if_else(n == 1, "Fusion", "Multi_Hit_Fusion"),
         # Required column for joining with MAF
         Variant_Type = "OTHER") %>%
  # Correct format for joining with MAF
  rename(Tumor_Sample_Barcode = Sample, Hugo_Symbol = Gene1) %>%
  # Create a new identifier from participant ID + sample_id
  ungroup() %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  rename(Tumor_Sample_Barcode = sample_id)

# Write to file
fusion_output <- file.path(output_dir, paste0(opt$filename_lead,
                                              "_fusions.tsv"))
readr::write_tsv(reformat_fusion, fusion_output)

#### CNV file preparation ------------------------------------------------------

cnv_df <- cnv_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  rename(Tumor_Sample_Barcode =  sample_id,
         Variant_Classification = status,
         Hugo_Symbol = gene_symbol) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification)

# Write to file
cnv_output <- file.path(output_dir, paste0(opt$filename_lead,
                                           "_cnv.tsv"))
readr::write_tsv(cnv_df, cnv_output)
