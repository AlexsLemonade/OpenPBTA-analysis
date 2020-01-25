# Format SV and CNV data into chromosomal breakpoint data
#
# C. Savonen for ALSF - CCDL
#
# 2020
#
# Code adapted from [svcnvplus](https://github.com/gonzolgarcia/svcnvplus).
#
# Option descriptions
# --cnv_seg: Relative file path (assuming from top directory of
#            'OpenPBTA-analysis') to CNV segment file to be used for breakpoint 
#            calculations.
# --sv: Relative file path (assuming from top directory of
#       'OpenPBTA-analysis') to SV file to be used for breakpoint calculations.
# --metadata: Relative file path (assuming from top directory of
#             'OpenPBTA-analysis') to MAF file to be analyzed. It is only needed
#             for identifying which samples are WGS or WXS.
# --ch.pct: A number between 0 and 1 that specifies the ratio of change needed to
#           consider a CNV copy number as changed.
# --output: Path to folder where you would like the output breakpoint calculations
#           files from this script to be stored
# --surveyed_wgs: File path that specifies the BED regions file that indicates
#                 the effectively surveyed regions of the genome for the WGS samples.",
# --surveyed_wxs: File path that specifies the BED regions file that indicates the
#                 effectively surveyed regions of the genome for the WXS samples.
#
# Command line example:
#
# Rscript analyses/chromosomal-instability/00-setup-breakpoint-data.R \
# --cnv_seg data/pbta-cnv-cnvkit.seg.gz \
# --sv data/pbta-sv-manta.tsv.gz \
# --metadata data/pbta-histologies.tsv \
# --surveyed_wgs  WGS_effectively_surveyed.bed \
# --surveyed_wxs data/WXS.hg38.100bp_padded.bed
#
################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = "--cnv_seg", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
    'OpenPBTA-analysis') to CNV segment file to be used for breakpoint 
    calculations.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--sv", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
    'OpenPBTA-analysis') to SV file to be used for breakpoint 
    calculations.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
    'OpenPBTA-analysis') to MAF file to be analyzed. It is only needed
    for identifying which samples are WGS or WXS. ",
    metavar = "character"
  ),
  make_option(
    opt_str = "--ch.pct", default = 0,
    help = "A number between 0 and 1 that specifies the ratio of change needed to
    consider a CNV copy number as changed.",
    metavar = "number"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output breakpoint calculations files from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--surveyed_wgs", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
    'OpenPBTA-analysis') that specifies the BED regions file that indicates the 
    effectively surveyed regions of the genome for the WGS samples.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--surveyed_wxs", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
    'OpenPBTA-analysis') that specifies the BED regions file that indicates the 
    effectively surveyed regions of the genome for the WXS samples.",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make everything relative to root path
opt$cnv_seg <- file.path(root_dir, opt$cnv_seg)
opt$sv <- file.path(root_dir, opt$sv)
opt$metadata <- file.path(root_dir, opt$metadata)
opt$surveyed_wgs <- file.path(root_dir, opt$surveyed_wgs)
opt$surveyed_wxs <- file.path(root_dir, opt$surveyed_wxs)

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$cnv_seg, opt$sv, opt$metadata, opt$surveyed_wgs, opt$surveyed_wxs
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
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

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


# Filter CNV data to only the changes that are larger than our cutoff `ch.pct`.
cnv_filtered_df <- cnv_df %>%
  dplyr::filter(changed == "change")

########################### Read in the SV data ################################
# Print progress message
message("Setting up SV breakpoints...")

sv_df <- data.table::fread(opt$sv, data.table = FALSE) %>%
  # Reformat the 23 and 24 chromosomes so they are X and Y and also factors
  dplyr::mutate(
    chrom = dplyr::recode(SV.chrom,
      `23` = "X", `24` = "Y"
    )
  )

#################### Format the data as chromosomes breaks #####################
# Only keep samples for which there are both SV and CNV data
common_samples <- dplyr::intersect(
  unique(cnv_df$ID),
  unique(sv_df$Kids.First.Biospecimen.ID.Tumor)
)

# TODO: after updating to consensus CNV data, evaluate whether sex chromosomes should still be removed.
# Make an CNV breaks data.frame.

# Make a CNV data.frame that has the breaks
cnv_breaks <- data.frame(
  samples = rep(cnv_filtered_df$ID, 2),
  chrom = rep(cnv_filtered_df$chrom, 2),
  coord = c(cnv_filtered_df$start, cnv_filtered_df$end),
  seg.mean = rep(cnv_filtered_df$seg.mean, 2),
  copy.num = rep(cnv_filtered_df$copy.num, 2),
  stringsAsFactors = FALSE
) %>%
  # Remove sex chromosomes
  dplyr::filter(
    !(chrom %in% c("X", "Y")),
    # Only keep samples that have both CNV and SV data
    samples %in% common_samples
  )

# Make an SV breaks data.frame.
sv_breaks <- data.frame(
  samples = rep(sv_df$Kids.First.Biospecimen.ID.Tumor, 2),
  chrom = rep(sv_df$chrom, 2),
  coord = c(sv_df$SV.start, sv_df$SV.end),
  svclass = rep(sv_df$SV.type, 2),
  stringsAsFactors = FALSE
) %>%
  # Remove sex chromosomes and NAs
  dplyr::filter(
    !(chrom %in% c("X", "Y", "M", NA)),
    # Only keep samples that have both CNV and SV data
    samples %in% common_samples
  )

############################## Create union of breaks ##########################
# Make an union of breaks data.frame.
union_of_breaks <- dplyr::bind_rows(sv_breaks, cnv_breaks) %>%
  dplyr::distinct(samples, chrom, coord, .keep_all = TRUE)

# Put all the breaks into a list.
breaks_list <- list(
  union_of_breaks = union_of_breaks,
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

# Set up the metadata
metadata <- readr::read_tsv(opt$metadata) %>%
  # Isolate metadata to only the samples that are in the datasets.
  dplyr::filter(Kids_First_Biospecimen_ID %in% common_samples) %>%
  # Keep the columns to only the experimental strategy and the biospecimen ID
  dplyr::select(Kids_First_Biospecimen_ID, experimental_strategy) %>%
  # For an easier time matching to our breaks data.frames, lets just rename this.
  dplyr::rename(samples = Kids_First_Biospecimen_ID)

# Calculate the breaks density for each data.frame
breaks_density_list <- lapply(breaks_list, function(breaks_df) {
  # Calculate the breaks density
  breaks_df %>%
    # Tack on the experimental strategy
    dplyr::inner_join(metadata) %>%
    # Recode using the BED range sizes
    dplyr::mutate(genome_size = dplyr::recode(experimental_strategy,
      "WGS" = wgs_size,
      "WXS" = wxs_size
    )) %>%
    dplyr::group_by(
      samples, experimental_strategy, genome_size
    ) %>%
    # Count number of mutations for that sample
    dplyr::summarize(breaks_count = dplyr::n()) %>%
    # Calculate breaks density
    dplyr::mutate(breaks_density = breaks_count / (genome_size / 1000000))
})

# Write the break densities each as their own files
purrr::imap(breaks_density_list, function(.x, name = .y) {
  # Write to TSV file
  readr::write_tsv(
    .x,
    file.path(opt$output, paste0(name, "_densities.tsv"))
  )
})
