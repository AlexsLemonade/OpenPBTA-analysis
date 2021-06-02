# J. Taroni for ALSF CCDL 2020
#
# Split up MAF files based on experimental strategy (WXS/WGS) for de novo
# signature extraction, per the sigfit vignette:
#
#  > For instance, if the cohort of samples includes both whole genomes and
#  > exomes, it is normally advisable to extract signatures from the
#  > whole-genome samples (provided that there are enough of them) and re-fit
#  > the resulting signatures to the whole-exome samples.
#
#  https://github.com/kgori/sigfit/blob/209776ee1d2193ad4b682b2e2472f848bd7c67a6/vignettes/sigfit_vignette.Rmd#L375
#
# Our cohort has both WXS and WGS samples.
#
# Everything in this script is harcoded since it is be highly specific to this
# project and its structure, naming, etc.

# Pipes
library(magrittr)

# Directory set up such that this can be run from anywhere within the project
# with hardcoding
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# Subdirectory of scratch for this module specifically, which we will create
# if it doesn't exist
scratch_dir <- file.path(root_dir, "scratch", "mutational-signatures")
if (!dir.exists(scratch_dir)) {
  dir.create(scratch_dir, recursive = TRUE)
}

# Read in the MAF file, subsetting to columns that are used by
# deconstructSigs::mut.to.sigs.input()
maf_file <- file.path(data_dir, "pbta-snv-consensus-mutation.maf.tsv.gz")
maf_df <- data.table::fread(maf_file, data.table = FALSE) %>%
  dplyr::select(Tumor_Sample_Barcode,
                Chromosome,
                Start_Position,
                Reference_Allele,
                Allele)

# Read in the metadata file, we only need the biospecimen ID which is the 
# identifier in the MAF-like file above
metadata_file <- file.path(data_dir, "pbta-histologies.tsv")
metadata_df <- readr::read_tsv(metadata_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                experimental_strategy)

# Join the experimental strategy to the relevant mutation columns
maf_df <- maf_df %>%
  dplyr::inner_join(metadata_df,
                    by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID"))

# Filter to WGS and write to file
wgs_file <- file.path(scratch_dir, "pbta-snv-consensus-wgs.tsv.gz")
maf_df %>%
  dplyr::filter(experimental_strategy == "WGS") %>%
  dplyr::select(-experimental_strategy) %>%
  readr::write_tsv(wgs_file)

# Filter to WXS
wxs_file <- file.path(scratch_dir, "pbta-snv-consensus-wxs.tsv.gz")
maf_df %>%
  dplyr::filter(experimental_strategy == "WXS") %>%
  dplyr::select(-experimental_strategy) %>%
  readr::write_tsv(wxs_file)
