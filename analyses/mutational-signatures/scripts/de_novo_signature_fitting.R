# Script by J. Taroni for ALSF CCDL
# Adapted from code written by Anna R. Poetsch and Candace L. Savonen
# TODO: document purpose, usage
# 

# Required for command line arguments
library(optparse)

# Required for signature fitting
library(deconstructSigs)
library(sigfit)

# Set up command line options
option_list <- list(
  make_option(c("--maf_file"),
              type = "character",
              help = "MAF file to be used for mutational signature extraction"),
  make_option(c("--nsignatures_floor"),
              type = "integer",
              default = 10,
              help = "The lower bound for the number of signatures to be extracted"),
  make_option(c("--nsignatures_ceiling"),
              type = "integer",
              default = 10,
              help = "The upper bound for the number of signatures to be extracted; should be >= --nsignature_floor"),
  make_option(c("--num_iterations"),
              type = "integer",
              default = 1000,
              help = "Number of iterations for sampling"),
  make_option(c("--seed"),
              type = "integer",
              default = 42,
              help = "Seed"),
  make_option(c("--output_file"),
              type = "character",
              help = "Full path to output file (.RDS)")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list = option_list))

# Warn if the ceiling for the number of signatures is less than the floor
if (opt$nsignatures_ceiling < opt$nsignatures_floor) {
  warning("Upper bound for signatures is less than lower bound!")
}

# The nsignatures argument to sigfit::extract_signatures can be a range of
# values, so we get the vector from the nsignatures_floor and
# nsignatures_ceiling arguments where if these numbers are the same this will be
# a vector of length 1
num_signatures <- opt$nsignatures_floor:opt$nsignatures_ceiling

# Read in the MAF file
maf <- data.table::fread(opt$maf_file, data.table = FALSE)

# Prep mutations file for signature extraction
sigs_input <- deconstructSigs::mut.to.sigs.input(
  mut.ref = maf,
  sample.id = "Tumor_Sample_Barcode",
  chr = "Chromosome",
  pos = "Start_Position",
  ref = "Reference_Allele",
  alt = "Allele",
  bsg = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

# Extract signatures
extracted_signatures <- sigfit::extract_signatures(
  counts = sigs_input,
  nsignatures = num_signatures,
  iter = opt$num_iterations,
  seed = opt$seed
)

# Write the extracted signatures list to (compressed) RDS file
readr::write_rds(extracted_signatures, path = opt$output_file, compress = "gz")
