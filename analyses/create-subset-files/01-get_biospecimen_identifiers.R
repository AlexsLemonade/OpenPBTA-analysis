# J. Taroni for CCDL 2019
# This script takes a directory of OpenPBTA files to subset and produces a list
# of biospecimen IDs, saved as an RDS file, to use to subset the files for
# use in continuous integration.
#
# This list will have the following features, where each element is a vector
# of biospecimen IDs to be extracted from a file:
#
#   - Some number of biospecimen IDs that correspond to participant IDs that
#     are represented across experimental strategies. The number of participant
#     IDs used to accomplish this is specified with the num_matched parameter.
#     Poly-A samples will be somewhat overrepresented vs. the portion of total
#     RNA-seq samples they make up.
#   - Some number of biospecimen IDs that correspond to participant IDs that
#     are *not* represented across strategies but are present in the file under
#     consideration. This number will be 10% of num_matched.
#   - We stratify based on `reported_gender`.
#   - We include (and hardcode) a set of biospecimen IDs for samples that have
#     TP53 and NF1 mutations that meet the criteria in the tp53_nf1_module and
#     are represented in the stranded RNA-seq dataset.
#     See 00-enrich-positive-examples for more information.
#
# EXAMPLE USAGE:
#
#   Rscript analyses/create-subset-files/01-get_biospecimen_identifiers.R \
#     --data_directory data/release-v5-20190924 \
#     --output_file analyses/create-subset-files/biospecimen_ids_for_subset.RDS \
#     --supported_string "pbta-snv|pbta-cnv|pbta-fusion|pbta-isoform|pbta-sv|pbta-gene" \
#     --num_matched 25 \
#     --seed 2019

#### Library and functions -----------------------------------------------------

library(readr)
library(optparse)

`%>%` <- dplyr::`%>%`

get_biospecimen_ids <- function(filename, id_mapping_df) {
  # Given a supported OpenPBTA file, return the participant IDs corresponding
  # to the biospecimen IDs contained within that file
  #
  # Args:
  #   filename: the full path to a supported OpenPBTA file
  #   id_mapping_df: the data.frame that contains mapping between biospecimen
  #                  IDs and participant IDs
  #
  # Returns:
  #   a vector of unique participant IDs

  message(paste("Reading in", filename, "..."))

  # where the biospecimen IDs come from in each file depends on the file
  # type -- that is why we need all of this logic
  if (grepl("pbta-snv", filename)) {
    # all SNV variant files keep the biospecimen identifiers in a column called
    # 'Tumor_Sample_Barcode'
    # if the files have consensus in the name, the first line of the file does
    # not contain MAF version information
    if (grepl("consensus|hotspots", filename)) {
      snv_file <- data.table::fread(filename, data.table = FALSE)
    } else {
      snv_file <- data.table::fread(filename,
                                    skip = 1,  # skip version string
                                    data.table = FALSE)
    }
    # both kinds (original, consensus)
    biospecimen_ids <- unique(snv_file$Tumor_Sample_Barcode)
  } else if (grepl("pbta-cnv", filename)) {
    # the two CNV files now have different structures
    cnv_file <- read_tsv(filename)
    if (grepl("controlfreec", filename)) {
      biospecimen_ids <- unique(cnv_file$Kids_First_Biospecimen_ID)
    } else {
      biospecimen_ids <- unique(cnv_file$ID)
    }
  } else if (grepl("consensus_seg_annotated", filename)) {
    annotated_cn_file <- read_tsv(filename)
    biospecimen_ids <- unique(annotated_cn_file$biospecimen_id)
  } else if (grepl("pbta-fusion", filename)) {
    fusion_file <- read_tsv(filename)
    # the biospecimen IDs in the filtered/prioritize fusion list included with
    # the download are in a column called 'Sample'
    if (grepl("putative-oncogenic|bysample", filename)) {
      biospecimen_ids <- unique(fusion_file$Sample)
    } else {
      # the original files contain the relevant IDs in a column 'tumor_id'
      biospecimen_ids <- unique(fusion_file$tumor_id)
    }
  } else if (grepl("pbta-sv", filename)) {
    # in a column 'Kids.First.Biospecimen.ID.Tumor'
    sv_file <- data.table::fread(filename, data.table = FALSE)
    biospecimen_ids <- unique(sv_file$Kids.First.Biospecimen.ID.Tumor)
  } else if (grepl(".rds", filename)) {
    # any column name that contains 'BS_' is a biospecimen ID
    expression_file <- read_rds(filename) %>%
      dplyr::select(dplyr::contains("BS_"))
    biospecimen_ids <- unique(colnames(expression_file))
  } else if (grepl("cnv_consensus", filename)) {
    cnv_consensus <- read_tsv(filename)
    biospecimen_ids <- unique(cnv_consensus$Biospecimen)
  } else {
    # error-handling
    stop("File type unrecognized by 'get_biospecimen_ids'")
  }

  # map from biospecimen ID to participant IDs and return unique participant IDs
  participant_ids <-  id_mapping_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% biospecimen_ids) %>%
    dplyr::pull(Kids_First_Participant_ID)
  return(unique(participant_ids))

}

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  make_option(
    c("-d", "--data_directory"),
    type = "character",
    default = NULL,
    help = "directory that contains data files to subset",
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "output RDS file"
  ),
  make_option(
    c("-r", "--supported_string"),
    type = "character",
    default = "pbta-snv|pbta-cnv|pbta-fusion|pbta-isoform|pbta-sv|pbta-gene|consensus_seg_annotated",
    help = "string for pattern matching used to subset to only supported files"
  ),
  make_option(
    c("-p", "--num_matched"),
    type = "integer",
    default = 25,
    help = "number of matched participants",
    metavar = "integer"
  ),
  make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2019,
    help = "seed integer",
    metavar = "integer"
  ),
  make_option(
    c("-l", "--local"),
    type = "integer",
    default = 0,
    help = "0 or 1; setting to 1 will skip the larger MAF files for local testing"
  )
)

# Read the arguments passed
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Handle options for whether or not this is running locally on someone's laptop
if (opt$local == 0) {
  running_locally <- FALSE
} else if (opt$local == 1) {
  running_locally <- TRUE
} else {
  stop("--local must be 0 or 1!")
}

# set up required arguments: input directory and output file
data_directory <- opt$data_directory
output_file <- opt$output_file

# string that will be used for pattern matching to select the files for
# subsetting
supported_files_string <- opt$supported_string

# get numbers of matched participants
num_matched_participants <- opt$num_matched
num_matched_polya <- ceiling(0.2 * num_matched_participants)
num_matched_stranded <- num_matched_participants - num_matched_polya

# extra non-matched samples -- this is a realistic scenario and will test for
# brittleness of code up for review
num_nonmatched <- ceiling(0.1 * num_matched_participants)

# set the seed
set.seed(opt$seed)

#### Samples we need to include to run tp53_nf1_score --------------------------
# For more information, see the 00-enrich-positive-examples notebook

tp53_dnaseq <- c("BS_3E5C1PN1", "BS_KHSYAB3J", "BS_ZR75EKKX")
tp53_stranded <- c("BS_TM9MH0RP", "BS_4B0BAVTX", "BS_FCAKKZ20")
nf1_dnaseq <- c("BS_RXSBJ929", "BS_P3PF53V8", "BS_7M7JNG00")
nf1_stranded <- c("BS_WJ9H4NZZ", "BS_Z8YZ6QGS", "BS_E0N3PTPN")

#### Get IDs -------------------------------------------------------------------

# list all files we are interested in subsetting and can support
files_to_subset <- list.files(data_directory,
                              pattern = supported_files_string,
                              full.names = TRUE)

# there are 6 RSEM files per library strategy, which we will assume contain the
# same samples
polya_rsem_files <- files_to_subset[grep("rsem.*polya", files_to_subset)]
stranded_rsem_files <- files_to_subset[grep("rsem.*stranded", files_to_subset)]

# we're going to remove 5 out of six of each set of files and use those samples
# for each of the six files later
kept_polya_rsem <- polya_rsem_files[1]
kept_stranded_rsem <- stranded_rsem_files[1]

# drop the unneeded files for each
files_to_subset <-
  files_to_subset[which(!(files_to_subset %in% c(polya_rsem_files[-1],
                                                 stranded_rsem_files[-1])))]

# drop recurrently fused genes by histology file -- it is small enough to
# include the entire thing
files_to_subset <-
  files_to_subset[-grep("fused-genes-byhistology", files_to_subset)]

# drop GISTIC zipped file from this list -- there are many files that are not
# currently documented
# we'll include the entire zipped folder
files_to_subset <-
  files_to_subset[-grep("gistic.zip", files_to_subset)]

# if testing this locally, drop the 2 larger of the 4 MAF files
if (running_locally) {
  files_to_subset <- files_to_subset[-grep("vardict|mutect2", files_to_subset)]
}

# get the participant ID to biospecimen ID
id_gender_df <- read_tsv(file.path(data_directory, "pbta-histologies.tsv"), guess_max = 10000) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID,
                reported_gender) %>%
  dplyr::distinct()

# drop reported gender
id_mapping_df <- id_gender_df %>%
  dplyr::select(-reported_gender)

# for each file, extract the participant ID list by first obtaining the
# biospecimen IDs and then mapping back to
participant_id_list <- purrr::map(files_to_subset,
                                  ~ get_biospecimen_ids(.x, id_mapping_df)) %>%
  purrr::set_names(files_to_subset)

# explicitly perform garbage collection here
gc(verbose = FALSE)

# split up information for poly-A and stranded expression data vs. all else
polya_participant_list <- purrr::keep(
  participant_id_list,
  grepl("polya",
        names(participant_id_list))
)

stranded_participant_list <- purrr::keep(
  participant_id_list,
  grepl("stranded",
        names(participant_id_list))
)

other_strategy_participant_list <- purrr::discard(
  participant_id_list,
  grepl("polya|stranded",
        names(participant_id_list))
)

# find samples that are common to all files
polya_in_all <- purrr::reduce(polya_participant_list, intersect)
stranded_in_all <- purrr::reduce(stranded_participant_list, intersect)
other_strategy_in_all <- purrr::reduce(other_strategy_participant_list,
                                       intersect)

# find RNA-seq samples that have matched samples in other strategies
polya_matched <- intersect(polya_in_all, other_strategy_in_all)
stranded_matched <- intersect(stranded_in_all, other_strategy_in_all)

# find identifiers for matched participants!
# first consider the poly-A samples
polya_for_subset <- sample(polya_matched, num_matched_polya)

# we need to sample the stranded participants keeping the reported gender in
# mind -- currently this is only for the sex prediction from RNA-seq data
id_gender_df <- id_gender_df %>%
  dplyr::filter(Kids_First_Participant_ID %in% stranded_matched)

# get the number of samples for each reported gender - we'll split 54% males
# which is what the 'matched' cohort is like
num_male <- ceiling(0.54 * num_matched_stranded)
num_female <- num_matched_stranded - num_male

stranded_for_subset <- c(
  id_gender_df %>%
    dplyr::filter(reported_gender == "Male") %>%
    dplyr::pull(Kids_First_Participant_ID) %>%
    unique() %>%
    sample(num_male),
  id_gender_df %>%
    dplyr::filter(reported_gender == "Female") %>%
    dplyr::pull(Kids_First_Participant_ID) %>%
    unique() %>%
    sample(num_female)
)

# intersect with the other strategies
matched_for_subset <- purrr::map(participant_id_list,
                                 ~ intersect(.x, c(polya_for_subset,
                                                   stranded_for_subset)))

# get a list of samples that are not in these matched lists to use to subset
nonmatched_for_subset <-
  purrr::map(participant_id_list,
             ~ setdiff(.x, c(polya_matched, stranded_matched))) %>%
  purrr::map(~ sample(.x, num_nonmatched))

# combine the matched and nonmatched lists of ids for subsetting
participant_ids_for_subset <- purrr::map2(matched_for_subset,
                                          nonmatched_for_subset,
                                          c)

# map back to biospecimen ids
biospecimen_ids_for_subset <- purrr::map(
  participant_ids_for_subset,
  function(x) {
    id_mapping_df %>%
      dplyr::filter(Kids_First_Participant_ID %in% x) %>%
      dplyr::pull(Kids_First_Biospecimen_ID)
  }
)

# for each stranded instance, add in biospecimen IDs for samples we know have a
# positive example of NF1 mutation and TP53 for tp53_nf1_score
stranded_index <- stringr::str_which(names(biospecimen_ids_for_subset),
                                     "stranded")
biospecimen_ids_for_subset <- biospecimen_ids_for_subset %>%
  purrr::modify_at(stranded_index, ~ append(.x, c(tp53_stranded, nf1_stranded)))

# for each pbta-snv instance, add in biospecimen IDs for samples we know have a
# positive example of NF1 mutation and TP53 for tp53_nf1_score
pbta_snv_index <- stringr::str_which(names(biospecimen_ids_for_subset),
                                     "pbta-snv")
biospecimen_ids_for_subset <- biospecimen_ids_for_subset %>%
  purrr::modify_at(pbta_snv_index, ~ append(.x, c(tp53_dnaseq, nf1_dnaseq)))

# now for the other RSEM files, we need to use the same identifiers as the
# same file we included
polya_rsem_ids <-
  biospecimen_ids_for_subset[[grep(kept_polya_rsem,
                                   names(biospecimen_ids_for_subset))]]
stranded_rsem_ids <-
  biospecimen_ids_for_subset[[grep(kept_stranded_rsem,
                                   names(biospecimen_ids_for_subset))]]

# create lists that contain the same identifiers
rest_polya_rsem <- lapply(polya_rsem_files[-1], function(x) polya_rsem_ids) %>%
  purrr::set_names(polya_rsem_files[-1])
rest_stranded_rsem <- lapply(stranded_rsem_files[-1],
                             function(x) stranded_rsem_ids) %>%
  purrr::set_names(stranded_rsem_files[-1])

# append the other RSEM elements to the list of all ids and write to file
biospecimen_ids_for_subset %>%
  append(rest_polya_rsem) %>%
  append(rest_stranded_rsem) %>%
  write_rds(output_file)
