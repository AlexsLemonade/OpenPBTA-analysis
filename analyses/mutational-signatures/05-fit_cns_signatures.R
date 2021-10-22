# Load libraries -------------------------------------------
library(optparse)
library(signature.tools.lib) # for getting the CNS signatures

# Set up command line options -------------------------------
option_list <- list(
  make_option(c("--abbreviated"),
              type = "integer", # should be integer or plays poorly with CI shell script parameter
              default = FALSE,
              action =  "store_true",
              help = "Run an abbreviated analysis with fewer iterations?"))
              

# Parse and check command line options ----------------------
opt <- parse_args(OptionParser(option_list = option_list))


# Full or abbreviated?
if (opt$abbreviated ==1) {
  n_iter <- 10 # CI time-saver
} else {
  n_iter <- 3000 # Full analysis
}

# Set up maf file and output file

# Path to root of project:
proj_root_path <- file.path( rprojroot::find_root(rprojroot::has_dir(".git")) )
analysis_path <- file.path(proj_root_path, "analyses", "mutational-signatures")

maf_file <- file.path(proj_root_path, "scratch", "mutational-signatures", "pbta-snv-consensus-wgs.tsv.gz")
fit_object_file <- file.path(analysis_path, "results", "fit_signatures_cns.rds")
fitted_exposures_file <- file.path(analysis_path, "results", "fitted_cns_signature_exposures.RDS")

# CNS signatures 
refsig_cns_matrix <- t( getOrganSignatures("CNS") )

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
fit <- sigfit::fit_signatures(counts = sigs_input, 
                              signatures = refsig_cns_matrix,
                              iter = n_iter, 
                              chains = 1, 
                              seed = 42, 
                              model = "poisson") # Highly similar performance to nmf and MUCH more efficient

# Save entire object, to be ignored
readr::write_rds(fit, fit_object_file, compress = "gz")


# Extract the information and save
fitted_exposures <- sigfit::retrieve_pars(fit, par = "exposures") 
readr::write_rds(fitted_exposures, fitted_exposures_file)

