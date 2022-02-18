# Load libraries -------------------------------------------
library(optparse)
library(signature.tools.lib) # contains the signal signatures
`%>%` <- dplyr::`%>%`


# Set up command line options -------------------------------
option_list <- list(
  make_option(c("--abbreviated"),
              type = "integer", # should be integer or plays poorly with CI shell script parameter
              default = FALSE,
              action =  "store_true",
              help = "Run an abbreviated analysis with fewer iterations for `sigfit`? This arg will be _ignored_ if --method is not `sigfit`.")
)


# Parse and check command line options ----------------------
opt <- parse_args(OptionParser(option_list = option_list))

# Full or abbreviated for MCMC?
if (opt$abbreviated == 1) {
  n_iter <- 10 # CI time-saver
} else {
  n_iter <- 3000 # Full analysis
}


# Set up directories, input/output files ----------------------------------

# Path to root of project:
proj_root_path <- file.path( rprojroot::find_root(rprojroot::has_dir(".git")) )
analysis_path <- file.path(proj_root_path, "analyses", "mutational-signatures")


# Path to input and output files:
maf_file <- file.path(proj_root_path, "data", "pbta-snv-consensus-mutation.maf.tsv.gz") # We consider _all_ WGS and WES mutations here
# Output files:
sigfit_fitted_file <- file.path(analysis_path, "results", "fitted_exposures_signal-cns_sigfit.rds")
decon_fitted_file  <- file.path(analysis_path, "results", "fitted_exposures_signal-cns_deconstructSigs.rds")


# Load CNS signatures 
refsig_cns_matrix <- t( getOrganSignatures("CNS") )

# Perform extraction -------------------------------
# Read in the MAF file
maf <- data.table::fread(maf_file, data.table = FALSE)

# Prep mutations file for signature fitting 
sigs_input <- deconstructSigs::mut.to.sigs.input(
  mut.ref = maf,
  sample.id = "Tumor_Sample_Barcode",
  chr = "Chromosome",
  pos = "Start_Position",
  ref = "Reference_Allele",
  alt = "Allele",
  bsg = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)


# Use sigfit to determine exposures to known CNS and save --------------------
fit1 <- sigfit::fit_signatures(counts = sigs_input, 
                               signatures = refsig_cns_matrix,
                               iter = n_iter, 
                               chains = 1, 
                               seed = 42, 
                               model = "poisson") # Highly similar performance to nmf and MUCH more efficient


fit1_exposures <- sigfit::retrieve_pars(fit1, par = "exposures") # extract exposures for saving
readr::write_rds(fit1_exposures, sigfit_fitted_file, compress = "gz")

# Use deconstructSigs to determine exposures to known CNS and save -------------------------
tumor_sample_ids <- maf %>%
  dplyr::filter(Tumor_Sample_Barcode %in% rownames(sigs_input)) %>%
  dplyr::distinct(Tumor_Sample_Barcode) %>%
  dplyr::pull(Tumor_Sample_Barcode)


fit2 <- lapply(tumor_sample_ids, function(sample_id) {
  #print(sample_id)
  # Determine the signatures contributing to the sample
  deconstructSigs::whichSignatures(
    tumor.ref = sigs_input,
    signatures.ref = as.data.frame(refsig_cns_matrix), # MUST BE DF!
    sample.id = sample_id,
    contexts.needed = TRUE
  )
})
readr::write_rds(fit2, decon_fitted_file, compress = "gz")


