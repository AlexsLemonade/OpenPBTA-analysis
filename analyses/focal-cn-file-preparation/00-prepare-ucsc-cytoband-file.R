# This script downloads and prepares the UCSC cytoband data to be added to the
# focal CN files that result from this module.
#
#
# Chante Bethell for CCDL 2020
#
# #### Example Usage
#
# This script is intended to be run via the command line.
# This example assumes it is being run from the root of the repository.
#
# Rscript --vanilla analyses/focal-cn-file-preparation/00-prepare-ucsc-cytoband-file.R

#### Set Up --------------------------------------------------------------------

# Install GenomicRanges
if (!("GenomicRanges" %in% installed.packages())) {
  BiocManager::install("GenomicRanges", update = FALSE)
}

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results directory
results_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#### Download and Wrangle UCSC data -------------------------------------------

# Read in UCSC cytoband data. The decision to implement the UCSC hg38 cytoband
# file was made based on a comparison done between the cytoband calls in the
# `org.Hs.eg.db` package and the calls in the UCSC file. We found that they
# disagreed in ~11,800 calls out of ~800,000 and the `UCSC file` contains more
# cytoband calls.
ucsc_cytoband <-
  data.table::fread(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
  )

# Make the UCSC cytoband data.frame a GRanges object
ucsc_cytoband_bed <- ucsc_cytoband %>%
  dplyr::select(chr = V1, start = V2, end = V3, cytoband = V4) %>%
  dplyr::mutate(cytoband = paste0(gsub("chr", "", chr), cytoband),
                chr = gsub("_.*","", chr)) %>%
  dplyr::filter(!(chr %in% c("chrUn", "chrM")))


readr::write_tsv(
  ucsc_cytoband_bed,
  file.path(results_dir, "ucsc_cytoband.bed")
)

#### Prepare consensus seg file -----------------------------------------------

consensus_with_status <-
  readr::read_tsv(file.path(root_dir, "scratch", "consensus_seg_with_status.tsv"))

consensus_with_status_bed <- consensus_with_status %>%
  dplyr::select(chr = chrom, start = loc.start, end = loc.end, status, Kids_First_Biospecimen_ID)

readr::write_tsv(
  consensus_with_status_bed,
  file.path(results_dir, "consensus_with_status.bed")
)
