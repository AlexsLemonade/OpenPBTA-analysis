# Load libraries -------------------------------------------
# Required for command line arguments
library(optparse)

# Required for signature extraction
library(deconstructSigs)
library(sigfit)
`%>%` <- dplyr::`%>%`


# Set up command line options -------------------------------
option_list <- list(
  make_option("--maf_file",
              type = "character",
              help = "MAF file to be used for mutational signature extraction"),    
  make_option("--seed",
              type = "integer",
              help = "The random seed to use for signature extraction"))


# Parse and check command line options ----------------------
opt <- parse_args(OptionParser(option_list = option_list))

# Check maf_file
maf_file <- opt$maf_file
if (!(file.exists(maf_file))) stop("MAF file does not exist at given path.")


 maf_file <- file.path("..", "..", "scratch", "mutational-signatures", "pbta-snv-consensus-wgs.tsv.gz")
# Check seed 
input_seed <- 1 #opt$seed

# CNS signatures matrix
path_to_refsig <- file.path("util", "refsig_cns_data", "refsig_cns_matrix.RDS")
refsig_cns_matrix <- readr::read_rds(path_to_refsig)



# Perform extraction -------------------------------
# Read in the MAF file
maf <- data.table::fread(maf_file, data.table = FALSE)

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


# Use sigfit to determine exposures to known CNS
fit <- fit_signatures(counts = sigs_input, 
                      signatures = refsig_cns_matrix,
                      iter = 3000, 
                      warmup = 1000, 
                      chains = 1, 
                      seed = input_seed, 
                      model = "poisson") # Highly similar performance to nmf and MUCH more efficient

# Extract the information and save here
fitted_signatures <- sigfit::retrieve_pars(fit, par = "signatures") 
fitted_exposures  <- sigfit::retrieve_pars(fit, par = "exposures") 

readr::write_rds(fitted_signatures, fitted_signatures_file)
readr::write_rds(fitted_exposures, fitted_exposures_file)


