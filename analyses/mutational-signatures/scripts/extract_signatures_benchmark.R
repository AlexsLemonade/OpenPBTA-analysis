# Written by Stephanie Spielman and adapted from code by Anna R. Poetsch, Candace L. Savonen, and J. Taroni
# Determine optimal k for de novo mutatational signature analysis as well as robustness to random seed.
#   Analysis is performed with a small number of iterations (1000) and tests k 2:8 (previous small runs suggest this range is appropraite) provided a given model (multinomial or poisson) and random seed.
#   Needs up to 20 GB RAM but is single-threaded
#   This script does NOT perform de novo signature extraction. It is the benchmark to determine parameters for actual inference.


# Load libraries -------------------------------------------
# Required for command line arguments
library(optparse)

# Required for signature extraction
library(deconstructSigs)
library(sigfit)



# Fixed benchmarking: We test k \in [2,8] with 1000 iterations each. More iterations is not sustainable for RAM.
num_signatures <- 2:8
num_iterations <- 1000

# Set up command line options -------------------------------
option_list <- list(
  make_option("--maf_file",
              type = "character",
              help = "MAF file to be used for mutational signature extraction"),    
  make_option("--model",
              type = "character",
              help = "The inference model for signature extraction, one of multinomial or poisson"),
  make_option("--seed",
              type = "integer",
              help = "The random seed to use for signature extraction"),
  make_option("--plot_dir",
              type = "character",
              help = "The output path to save goodness-of-fit plots to."))
                
                            
# Parse and check command line options ----------------------
opt <- parse_args(OptionParser(option_list = option_list))

# Check maf_file
maf_file <- opt$maf_file
if (!(file.exists(maf_file))) stop("MAF file does not exist at given path.")

# Check model
model <- opt$model
if (!(model %in% c("poisson", "multinomial"))) stop("Model must be one of: multinomial poisson")

# Check seed
seed <- opt$seed

# Make output_path if doesn't exist
plot_dir <- opt$plot_dir
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Use opts to define final PDF goodness-of-fit file
output_gof_plot <- file.path(plot_dir,
                             paste0("gof_seed_",seed, "_model_", model, ".pdf"))
print(output_gof_plot)


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

# Extract signatures
extracted_signatures <- sigfit::extract_signatures(
  counts      = sigs_input,
  nsignatures = num_signatures,
  iter        = num_iterations,
  seed        = seed,
  model       = model
)

# Save goodness-of-fit plot, but we do not save the extracted signature output since they are massive.
pdf(output_gof_plot, width = 7, height = 7)
sigfit::plot_gof(extracted_signatures)
dev.off()



































