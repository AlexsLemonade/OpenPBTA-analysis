# Format SV and CNV data into chromosomal breakpoint data
#
# C. Savonen for ALSF - CCDL
#
# 2020
#
# Code adapted from [svcnvplus](https://github.com/gonzolgarcia/svcnvplus).
#
# Option descriptions
# --cnv_seg: File path to CNV segment file to be used for breakpoint
#            calculations.
# --sv: File path to SV file to be used for breakpoint calculations.
# --metadata: File path to metadata file to be analyzed. It is only needed
#             for identifying which samples are WGS or WXS.
# --ch.pct: A number between 0 and 1 that specifies the ratio of change needed to
#           consider a CNV copy number as changed.
# --output: Path to folder where you would like the output breakpoint calculations
#           files from this script to be stored
# --surveyed_wgs: File path that specifies the BED regions file that indicates
#                 the effectively surveyed regions of the genome for the WGS samples.
# --surveyed_wxs: File path that specifies the BED regions file that indicates the
#                 effectively surveyed regions of the genome for the WXS samples.
# --gap: An integer that indicates how many base pairs away a CNV and SV and
#        still be considered the same. Will be passed to maxgap argument in
#        GenomicRanges::findOverlaps. Default is 0.
# --drop_sex: If TRUE, will drop the sex chromosomes. Default is FALSE
#
# Command line example:
#
# Rscript analyses/chromosomal-instability/00-setup-breakpoint-data.R \
# --cnv_seg data/pbta-cnv-cnvkit.seg.gz \
# --sv data/pbta-sv-manta.tsv.gz \
# --metadata data/pbta-histologies.tsv \
# --surveyed_wgs  scratch/WGS_effectively_surveyed.bed \
# --surveyed_wxs data/WXS.hg38.100bp_padded.bed
#
################################ Initial Set Up ################################
# Establish base dir, this is required for correctly sourcing the functions
# below
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# We need the `make_granges` function from here
source(file.path(root_dir,
                 "analyses",
                 "chromosomal-instability",
                 "util",
                 "chr-break-calculate.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)


############################ Intersect function ################################
intersect_cnv_sv <- function(sample_id, sv_breaks, cnv_breaks, gap = opt$gap) {
  # For a given sample's data in CNV and SV chromosomal breaks data.frame,
  # intersect the two (based on the maxgap allowed to consider two breaks
  # identical) and return a data.frame with the intersection breaks.
  #
  # Args:
  #   sample_id: The sample_id to be looked up in the samples_col
  #   sv/cnv_breaks: for a data.frame with chromosomal coordinates and sample IDs
  #                  for their respective breaks
  #   sample_id: a character string that designates which sample's data needs to be
  #              extracted and intersected between the two data.frames (CNV and SV)
  #   gap : The max number of bases between a CNV and SV break for them to be
  #         considered the same.
  #
  # Returns:
  # A chromosomal breaks data.frame that contains the intersection of CNV and SV
  # chromosomal break data.
  #
  # Make into GenomicRanges objects
  sv_ranges <- make_granges(sv_breaks,
                            sample_id = sample_id,
                            start_col = "coord",
                            end_col = "coord")
  cnv_ranges <- make_granges(cnv_breaks,
                             sample_id = sample_id,
                             start_col = "coord",
                             end_col = "coord")
  # Find overlaps
  intersection_df <- IRanges::mergeByOverlaps(
    sv_ranges,
    cnv_ranges,
    maxgap = opt$gap) %>%
    # Coerce to data.frame
    as.data.frame() %>%
    # Remove these, they are redundant
    dplyr::select(-dplyr::contains("mcols"))

  return(intersection_df)
}

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = "--cnv_seg", type = "character", default = "none",
    help = "File path to CNV segment file to be used for breakpoint calculations.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--sv", type = "character", default = "none",
    help = "File path to SV file to be used for breakpoint calculations.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", default = "none",
    help = "File path to MAF file to be analyzed. It is only needed for identifying which samples are WGS or WXS.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--ch.pct", default = 0,
    help = "A number between 0 and 1 that specifies the ratio of change needed to consider a CNV copy number as changed.",
    metavar = "number"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL,
    help = "Path to folder where you would like the output breakpoint calculations files from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--surveyed_wgs", type = "character", default = "none",
    help = "File path that specifies the BED regions file that indicates the effectively surveyed regions of the genome for the WGS samples.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--surveyed_wxs", type = "character", default = "none",
    help = "File path that specifies the BED regions file that indicates the effectively surveyed regions of the genome for the WXS samples.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--uncalled_samples", type = "character", default = "none",
    help = "File path that specifies the regions that were uncalled in CNV analysis.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--gap", default = 0,
    help = "An integer that indicates how many base pairs away a CNV and SV breakpoint can be and still be considered the same. Will be passed to maxgap argument in GenomicRanges::findOverlaps.",
    metavar = "number"
  ),
  make_option(
    opt_str = "--drop_sex", action = "store_true",
    default = FALSE, help = "If TRUE, will drop the sex chromosomes. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$cnv_seg, opt$sv, opt$metadata, opt$surveyed_wgs, opt$surveyed_wxs, opt$uncalled
)

# Get list of which files were found
files_found <- file.exists(needed_files)

# Report error if any of them aren't found
if (!all(files_found)) {
  stop(paste("\n Could not find needed file(s):",
    needed_files[which(!files_found)],
    "Check your options and set up.",
    sep = "\n"
  ))
}

############################## Set Up Output ###################################

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

################################ Set up CNV data ###############################
# Print progress message
message("Setting up CNV breakpoints...")

# Read in the segment copy number data
cnv_df <- data.table::fread(opt$cnv_seg, data.table = FALSE) %>%
  dplyr::mutate(
    samples = as.factor(ID),
    # Calculate width
    width = loc.end - loc.start,
    # Make a binary variable that labels whether or not a sequence
    # is considered changed by the threshold set
    # TODO: after updating to Consensus CNV data, see if sex chromosomes still
    # need to be removed and whether these thresholds should be changed.
    changed = as.factor(dplyr::case_when(
      abs(seg.mean) > log2(1 + opt$ch.pct) ~ "change",
      TRUE ~ "no_change"
    ))
  ) %>%
  # Reformat the chromosome variable to drop the "chr"
  dplyr::mutate(chrom = factor(gsub("chr", "", chrom),
    levels = c(1:22, "X", "Y")
  )) %>%
  # Changing these so they end up matching the SV data
  dplyr::rename(start = loc.start, end = loc.end)

# Obtain full list of samples
cnv_samples <- unique(as.character(cnv_df$samples))

# Filter CNV data to only the changes that are larger than our cutoff `ch.pct`.
cnv_filtered_df <- cnv_df %>%
  dplyr::filter(changed == "change")

########################### Read in the SV data ################################
# Print progress message
message("Setting up SV breakpoints...")

sv_df <- data.table::fread(opt$sv, data.table = FALSE) %>%
  # Filter for PASS variants only
  dplyr::filter(FILTER == "PASS") %>%
  # Reformat the 23 and 24 chromosomes so they are X and Y and also factors
  dplyr::mutate(
    chrom = dplyr::recode(SV.chrom,
      `23` = "X", `24` = "Y"
    )
  )

# Obtain full list of samples
sv_samples <- unique(sv_df$Kids.First.Biospecimen.ID.Tumor)

####################### Drop Sex Chr if option is on ###########################
if (opt$drop_sex){
  sv_df <- sv_df %>%
    dplyr::filter(!(chrom %in% c("X", "Y", "M")))

  cnv_df <- cnv_df %>%
    dplyr::filter(!(chrom %in% c("X", "Y", "M")))
}
#################### Format the data as chromosomes breaks #####################
# Make a CNV data.frame that has the breaks
cnv_breaks <- data.frame(
  samples = rep(cnv_filtered_df$ID, 2),
  chrom = rep(cnv_filtered_df$chrom, 2),
  coord = c(cnv_filtered_df$start, cnv_filtered_df$end),
  seg.mean = rep(cnv_filtered_df$seg.mean, 2),
  copy.num = rep(cnv_filtered_df$copy.num, 2),
  stringsAsFactors = FALSE
) %>%
  # Remove NAs
  dplyr::filter(!is.na(chrom))

# Make an SV breaks data.frame.
sv_breaks <- data.frame(
  samples = rep(sv_df$Kids.First.Biospecimen.ID.Tumor, 2),
  chrom = rep(sv_df$chrom, 2),
  coord = c(sv_df$SV.start, sv_df$SV.end),
  svclass = rep(sv_df$SV.type, 2),
  stringsAsFactors = FALSE
) %>%
  # Remove NAs
  dplyr::filter(!is.na(chrom))

######################### Create intersection of breaks ########################
# Get list of samples for which there are both SV and CNV data
common_samples <- dplyr::intersect(
  unique(cnv_breaks$samples),
  unique(sv_breaks$samples)
)

# Make an intersection of breaks data.frame.
intersection_of_breaks <- lapply(common_samples,
                                 intersect_cnv_sv, # Special intersect function
                                 sv_breaks = sv_breaks,
                                 cnv_breaks = cnv_breaks,
                                 gap = opt$gap) # The maxgap allowed

# Bring along sample names
names(intersection_of_breaks) <- common_samples

# Collapse into a data.frame
intersection_of_breaks <- dplyr::bind_rows(intersection_of_breaks,
                                           .id = "samples") %>%
  # Only need one chromosome
  dplyr::select(chrom = sv_ranges.seqnames, dplyr::everything(), -cnv_ranges.seqnames)

# Put all the breaks into a list.
breaks_list <- list(
  intersection_of_breaks = intersection_of_breaks,
  cnv_breaks = cnv_breaks,
  sv_breaks = sv_breaks
)

# Save to an RDS file
readr::write_rds(breaks_list, file.path(opt$output, "breaks_lists.RDS"))

################# Calculate the breaks density using BED ranges ################
# Set up the BED ranges for the denominator of break density
# Read in the BED files we need
bed_wgs <- readr::read_tsv(opt$surveyed_wgs, col_names = FALSE)
bed_wxs <- readr::read_tsv(opt$surveyed_wxs, col_names = FALSE)

# Sum up genome sizes
wgs_size <- sum(bed_wgs[, 3] - bed_wgs[, 2])
wxs_size <- sum(bed_wxs[, 3] - bed_wxs[, 2])

# Don't want integers per se
wgs_size <- as.numeric(wgs_size)
wxs_size <- as.numeric(wxs_size)

# Read in uncalled samples from consensus module
uncalled_samples <- readr::read_tsv(opt$uncalled)$sample

# Get vector of all samples
all_samples <- unique(c(cnv_samples, sv_samples, uncalled_samples))


# Set up the metadata
metadata <- readr::read_tsv(opt$metadata, guess_max = 10000) %>%
  # Isolate metadata to only the samples that are in the datasets.
  dplyr::filter(Kids_First_Biospecimen_ID %in% all_samples) %>%
  # Keep the columns to only the experimental strategy and the biospecimen ID
  dplyr::select(Kids_First_Biospecimen_ID, experimental_strategy) %>%
  # For an easier time matching to our breaks data.frames, lets just rename this.
  dplyr::rename(samples = Kids_First_Biospecimen_ID) %>%
  # samples not in cnv_samples were were noisy or otherwise bad data
  dplyr::mutate(surveyed = samples %in% cnv_samples)


# Calculate the breaks density for each data.frame
breaks_density_list <- lapply(breaks_list, function(breaks_df) {
  # Calculate the breaks density
  breaks_df %>%
    # Tack on the experimental strategy and samples with no counts
    dplyr::full_join(metadata) %>%
    # Recode using the BED range sizes
    dplyr::mutate(genome_size = dplyr::recode(experimental_strategy,
      "WGS" = wgs_size,
      "WXS" = wxs_size
    )) %>%
    dplyr::group_by(
      samples, experimental_strategy, genome_size, surveyed
    ) %>%
    # Count number of mutations for that sample but find out if it is NA
    dplyr::summarize(is_na = any(is.na(chrom)),
                     breaks_count = dplyr::n()) %>%
    # Calculate breaks density, but put NA for breaks_count if the sample was
    # dropped from CNV consensus analysis
    dplyr::mutate(breaks_count = dplyr::case_when(
      !is_na ~ as.numeric(breaks_count),
      is_na & surveyed ~ as.numeric(0),
      TRUE ~ as.numeric(NA)
      ),
    breaks_density = breaks_count / (genome_size / 1000000)) %>%
    # Drop the is_na column, we only needed if for recoding
    dplyr::select(-is_na, -surveyed)
})

# Write the break densities each as their own files
purrr::imap(breaks_density_list, function(.x, name = .y) {
  # Write to TSV file
  readr::write_tsv(
    .x,
    file.path(opt$output, paste0(name, "_densities.tsv"))
  )
})
