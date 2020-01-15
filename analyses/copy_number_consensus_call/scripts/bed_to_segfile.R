
# Base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "copy_number_consensus_call")

library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Parse options

option_list <- list(
  make_option(
    c("-i", "--cnv_file"),
    type = "character",
    default = NULL,
    help = "path to the cnv consensus file, relative to analysis root",
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "path to output file, relative to analysis root"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# read the cnv output table
cnvs <- readr::read_tsv(file.path(analysis_dir, opts$cnv_file))

# Columns *_CNV have the format START:END:COPY_NUMBER:SEG.MEAN,START:END:COPY_NUMBER:SEG.MEAN
# This function converts one into a data frane, 
# which we can then make a column within the cnv data frame
segstrings_to_df <- function(segment_string){
  segment_vector <- stringr::str_split(segment_string, ",")[[1]]
  seg_matrix <- stringr::str_split(segment_vector, ":", simplify = TRUE)
  # account for rows with no calls
  if(dim(seg_matrix)[2] != 4){seg_matrix <- matrix(ncol = 4)}
  
  colnames(seg_matrix) <- c("start", "end", "copies", "segmean")
  suppressWarnings(class(seg_matrix) <- "numeric")
  seg_df <- as.data.frame(seg_matrix) %>%
    dplyr::mutate(length = abs(end - start))
  return(seg_df)
}

# Functions to calculate weighted means and medians 
# from data frames created with segstrings_to_df()

seg_wmean <- function(seg_df){
  wmean <- weighted.mean(seg_df$segmean, seg_df$length)
}

seg_wmedian <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$segmean, seg_df$length, 
                                      interpolate = FALSE)
}

seg_wmedian2 <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$segmean, seg_df$length, 
                                      interpolate = TRUE)
}

copies_wmedian <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$copies, seg_df$length, 
                                      interpolate = FALSE)
}


# Calculate summary stats from merged CNV calls. \
cnvs <- cnvs %>%
  dplyr::mutate(cnvkit_df = purrr::map(cnvkit_CNVs, segstrings_to_df),
                freec_df = purrr::map(freec_CNVs, segstrings_to_df),
                segmean_wmean = purrr::map_dbl(cnvkit_df, seg_wmean),
                segmean_wmedian1 = purrr::map_dbl(cnvkit_df, seg_wmedian),
                segmean_wmedian2 = purrr::map_dbl(cnvkit_df, seg_wmedian2),
                cnvkit_cn = purrr::map_dbl(cnvkit_df, copies_wmedian),
                freec_cn = purrr::map_dbl(freec_df, copies_wmedian),
                copynum = ifelse(is.finite(cnvkit_cn), # use cnvkit if available
                                 cnvkit_cn, freec_cn),
                num.mark = NA)


# .seg output format from CNVkit is: 
# ID	chrom	loc.start	loc.end	num.mark	seg.mean	copy.num

out_table <- cnvs %>%
  dplyr::select(ID = Biospecimen,
                chrom = chrom,
                loc.start = start, 
                loc.end = end,
                num.mark = num.mark,
                seg.mean = segmean_wmean)

readr::write_tsv(out_table, file.path(analysis_dir, opts$output_file))