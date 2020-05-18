# J. Taroni for ALSF CCDL 2019
# This script processes MAF, focal CN (from focal-cn-file-preparation),
# standardized fusion files and prepares them for oncoprint plotting.
#
# NOTES:
#   * The `Tumor_Sample_Barcode` will now corresponds to the `sample_id` column
#     in the histologies file
#   * We remove ambiguous `sample_id` -- i.e., where there are more than two
#     tumor biospecimens that map to the same sample id
#   * Filtering via an independent specimen file is optional, but highly
#     recommended
#
# EXAMPLE USAGE:
#
# Rscript --vanilla 00-map-to-participant.R \
#   --maf_file snv-consensus_11122019/consensus_mutation.maf.tsv \
#   --cnv_file ../focal-cn-file-preparation/results/controlfreec_annotated_cn_autosomes.tsv.gz \
#   --fusion_file ../../scratch/arriba.tsv \
#   --metadata_file ../../data/pbta-histologies.tsv \
#   --output_directory ../../scratch/oncoprint_files \
#   --filename_lead "all_participants_primary_only" \
#   --independent_specimens ../../data/independent-specimens.wgswxs.primary.tsv

library(dplyr)

#### Command line options ------------------------------------------------------

option_list <- list(
  optparse::make_option(
    c("--maf_file"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains SNV information",
  ),
  optparse::make_option(
    c("--cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains CNV information"
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
    help = "character string used at the beginning of output file names"
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

#### Directory and output file set up ------------------------------------------

output_dir <- opt$output_directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# All of the output files will use the filename lead supplied as an argument
# followed by what kind of data they contain
maf_output <- file.path(output_dir, paste0(opt$filename_lead, "_maf.tsv"))
fusion_output <- file.path(output_dir, paste0(opt$filename_lead,
                                              "_fusions.tsv"))
cnv_output <- file.path(output_dir, paste0(opt$filename_lead, "_cnv.tsv"))

#### Read in data --------------------------------------------------------------

histologies_df <- readr::read_tsv(opt$metadata_file)
maf_df <- readr::read_tsv(opt$maf_file)
cnv_df <- readr::read_tsv(opt$cnv_file)
fusion_df <- readr::read_tsv(opt$fusion_file)

#### Get rid of ambiguous and non-tumor samples --------------------------------

# An ambiguous sample_id will have more than 2 rows associated with it in the
# histologies file when looking at tumor samples -- that means we won't be able
# to determine when an WGS/WXS assay maps to an RNA-seq assay for the purpose of
# the oncoprint plot
ambiguous_sample_ids <- histologies_df %>%
  filter(sample_type == "Tumor",
         composition == "Solid Tissue") %>%
  group_by(sample_id) %>%
  tally() %>%
  filter(n > 2) %>%
  pull(sample_id)

ambiguous_biospecimens <- histologies_df %>%
  filter(sample_id %in% ambiguous_sample_ids) %>%
  pull(Kids_First_Biospecimen_ID)

# we're only going to look at tumor samples in the oncoprint plot
not_tumor_biospecimens <- histologies_df %>%
  filter(sample_type != "Tumor",
         composition != "Solid Tissue") %>%
  pull(Kids_First_Biospecimen_ID)

biospecimens_to_remove <- unique(c(ambiguous_biospecimens,
                                   not_tumor_biospecimens))

# Filter the files!
maf_df <- maf_df %>%
  dplyr::filter(!(Tumor_Sample_Barcode %in% biospecimens_to_remove))
cnv_df <- cnv_df %>%
  dplyr::filter(!(biospecimen_id %in% biospecimens_to_remove))
fusion_df <- fusion_df %>%
  dplyr::filter(!(Sample %in% biospecimens_to_remove))

#### Filter to independent specimens (optional) --------------------------------

if (!is.null(opt$independent_specimens)) {

  ind_biospecimen <- readr::read_tsv(opt$independent_specimens) %>%
    pull(Kids_First_Biospecimen_ID)

  # filter the genome data, e.g., SNV and CNV data, to only include biospecimen
  # identifiers in the independent file
  maf_df <- maf_df %>%
    filter(Tumor_Sample_Barcode %in% ind_biospecimen)
  cnv_df <- cnv_df %>%
    filter(biospecimen_id %in% ind_biospecimen)

  # for the RNA-seq samples, we need to map from the sample identifier
  # associated with the independent specimen and back to a biospecimen ID
  ind_sample_id <- histologies_df %>%
    filter(Kids_First_Biospecimen_ID %in% ind_biospecimen) %>%
    pull(sample_id)

  # get the corresponding biospecimen ID to be used
  rnaseq_ind <- histologies_df %>%
    filter(sample_id %in% ind_sample_id,
           experimental_strategy == "RNA-Seq") %>%
    pull(Kids_First_Biospecimen_ID)

  # finally filter the fusions
  fusion_df <- fusion_df %>%
    filter(Sample %in% rnaseq_ind)

}

#### MAF file preparation ------------------------------------------------------

message("Preparing MAF file...")

# join the sample_id information to the MAF file and then set as the tumor
# sample barcode
maf_df <- maf_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  rename(Tumor_Sample_Barcode = sample_id)

# Write MAF to file
readr::write_tsv(maf_df, maf_output)

#### Fusion file preparation ---------------------------------------------------
# TODO: Once the consensus calls of the fusion data are obtained, this section
# will need to be adapted to the format of the fusion input file. For example,
# the way we separate the genes out of `FusionName` may need to be adapted.

message("Preparing fusion file...")
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
  select(-n) %>%
  ungroup() %>%
  # Correct format for joining with MAF
  rename(Tumor_Sample_Barcode = Sample, Hugo_Symbol = Gene1) %>%
  # Create a new identifier from participant ID + sample_id
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  mutate(Tumor_Sample_Barcode = sample_id)

# Write to file
readr::write_tsv(reformat_fusion, fusion_output)

#### CNV file preparation ------------------------------------------------------

message("Preparing CNV file...")
cnv_df <- cnv_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  mutate(Tumor_Sample_Barcode =  sample_id) %>%
  rename(Variant_Classification = status,
         Hugo_Symbol = gene_symbol) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification)

# Write to file
readr::write_tsv(cnv_df, cnv_output)
