# Script by J. Taroni for ALSF CCDL
# Adapted from code written by Anna R. Poetsch and Candace L. Savonen
# 
# Given a MAF file, extract de novo mutational signatures for a range of k as
# specified by the nsignatures_floor and nsignatures_ceiling arguments using 
# the sigfit package.
# 
# This script is essentially a wrapper for sigfit::extract_signatures().
# 
# Additional arguments control the number of iterations, the seed, the output
# RDS file and, optionally, the full path for a PDF of the goodness-of-fit plot
# from sigfit::plot_gof().
# 
# COMMAND LINE USAGE:
# 
### To run for a range of values to first identify the number of signatures to
### use, using a low number of iterations
#
# Rscript --vanilla \
#   scripts/de_novo_signature_extraction.R \
#   --maf_file ../../data/pbta-snv-consensus-mutation.maf.tsv.gz \
#   --nsignatures_floor 5 \
#   --nsignatures_ceiling 15 \
#   --num_iterations 1000 \
#   --seed 42 \
#   --output_file results/denovo_sigfit_signatures.RDS
#   
### To run once the number of signatures has been selected, in this example 10,
### with an appropriate number of iterations
#
# Rscript --vanilla \
#   scripts/de_novo_signature_extraction.R \
#   --maf_file ../../data/pbta-snv-consensus-mutation.maf.tsv.gz \
#   --nsignatures_floor 10 \
#   --nsignatures_ceiling 10 \
#   --num_iterations 10000 \
#   --seed 42 \
#   --output_file results/nsig10_signatures.RDS

#### ---------------------------------------------------------------------------

# Required for command line arguments
library(optparse)

# Required for signature extraction
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
  make_option("--model",
              type = "character",
              default = "multinomial",
              help = "The inference model for signature extraction, one of multinomial or poisson"),
  make_option(c("--seed"),
              type = "integer",
              default = 42,
              help = "Seed"),
  make_option(c("--output_file"),
              type = "character",
              default = NULL,
              help = "Full path to output file (.RDS)"),
  make_option(c("--plot_output"),
              type = "character",
              default = NULL,
              help = "Full path to goodness-of-fit plot (PDF)")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list = option_list))

# Check maf_file
if (!(file.exists(opt$maf_file))) stop("MAF file does not exist at given path.")

# Check model
if (!(opt$model %in% c("poisson", "multinomial"))) stop("Model must be one of: multinomial poisson")

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
  model = opt$model,
  iter = opt$num_iterations,
  seed = opt$seed
)

# If user specifies a plot output, save the goodness-of-fit plot at that 
# location
if (!is.null(opt$plot_output)) {
  png(opt$plot_output, 640, 480)
  sigfit::plot_gof(extracted_signatures)
  dev.off()
}

# If user specifies a file output, save the extracted signatures list
if (!is.null(opt$output_file)) {
  # Write the extracted signatures list to (compressed) RDS file
  readr::write_rds(extracted_signatures, path = opt$output_file, compress = "gz")
}
