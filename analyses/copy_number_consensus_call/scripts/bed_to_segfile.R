# bed_to_segfile.R
#
# Josh Shapiro for CCDL 2020
# 
# Purpose: Convert the bed file output from the CNV consensus workflow to a seg file
#
# Option descriptions
# -i, --cnv_file : path to the cnv consensus file 
# -o, --output_file : path for output file
# --segmean-method : method for combining seg.mean values. Default `weight_mean` for the length-weighted mean
#
# example invocation:
# Rscript scripts/bed_to_segfile.R \
#   -i results/cnv_consensus.tsv \
#   -o results/pbta-cnv-consensus.seg



# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Parse command line options
library(optparse)
option_list <- list(
  make_option(
    c("-i", "--cnv_file"),
    type = "character",
    default = NULL,
    help = "path to the cnv consensus file",
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "path to output file"
  ),
  make_option(
    c("--segmean_method"),
    type = "character",
    
    default = "weight_mean",
    help = "Method for combining segment mean values.
                One of ['weight_mean', 'weight_median', 'median_interpolate'],
                corresponding to weighted mean, weighted median, or weighted  
                median with interpolation. Default is weighted mean."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Check options:
if (! (opts$segmean_method %in% c('weight_mean', 'weight_median', 'median_interpolate'))){
    stop("--segmean_method must be one of ['weight_mean', 'weight_median', 'median_interpolate']")
}

# Read the cnv output table
cnvs <- readr::read_tsv(opts$cnv_file)

# Columns *_CNV have the format START:END:COPY_NUMBER:SEG.MEAN,START:END:COPY_NUMBER:SEG.MEAN
# This function converts ones such string into a data frane, 
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

## Functions to calculate weighted means and medians 
## from data frames created with segstrings_to_df()

# weighted mean segmean
seg_wmean <- function(seg_df){
  wmean <- weighted.mean(seg_df$segmean, seg_df$length)
}

# weighted median segmean
seg_wmedian <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$segmean, seg_df$length, 
                                      interpolate = FALSE)
}
## weighted median with interpolation (probably not used)
seg_wmedian_interpolate <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$segmean, seg_df$length, 
                                      interpolate = TRUE)
}

# weighted median copy number
copies_wmedian <- function(seg_df){
  wmed <- matrixStats::weightedMedian(seg_df$copies, seg_df$length, 
                                      interpolate = FALSE)
}

segmean_function <- if (opts$segmean_method == "weight_mean"){
    seg_wmean
  } else if (opts$segmean_method == "weight_median"){
    seg_wmedian
  } else if (opts$segmean_method == "median_interpolate"){
    seg_wmedian_interpolate
  }

# Calculate summary stats from merged CNV calls. \
cnvs <- cnvs %>%
  dplyr::mutate(cnvkit_df = purrr::map(cnvkit_CNVs, segstrings_to_df),
                freec_df = purrr::map(freec_CNVs, segstrings_to_df),
                segmean = purrr::map_dbl(cnvkit_df, segmean_function),
                cnvkit_cn = purrr::map_dbl(cnvkit_df, copies_wmedian),
                freec_cn = purrr::map_dbl(freec_df, copies_wmedian),
                copynum = ifelse(is.finite(cnvkit_cn), # use cnvkit if available
                                 cnvkit_cn, freec_cn), #otherwise use freec
                num.mark = NA)


# .seg output format from CNVkit is: 
# ID	chrom	loc.start	loc.end	num.mark	seg.mean	copy.num

out_table <- cnvs %>%
  dplyr::select(ID = Biospecimen,
                chrom = chrom,
                loc.start = start, 
                loc.end = end,
                num.mark = num.mark,
                seg.mean = segmean,
                copy.num = copynum)

readr::write_tsv(out_table, opts$output_file)
