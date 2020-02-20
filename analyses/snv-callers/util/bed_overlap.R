# Functions for calculating tumor mutation burden
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

bed_to_granges <- function(bed_file) {
  # For a given breaks data.frame make a GenomicRanges object from it. Can
  # filter to a single samples' data.
  #
  # Args:
  #   bed_file:
  #
  # Returns:
  # A BED file is read in and made into a GenomicRanges object.

  # If no bed_file is provided, stop
  if (file.exists(bed_file)) {
    stop("No BED file has been found at this file path. ")
  }

  bed_df <- readr::read_tsv(bed_file, col_names = c("chr", "start", "end"),
  col_types = list(col_character(), col_double(), col_double())) %>%
    dplyr::filter(!is.na(chr)) %>%
    dplyr::mutate(chr = stringr::word(chr, sep = "_", 1))

  # Stop if the sample doesn't exist in the data.frame
  if (nrow(bed_df) == 0) {
    stop("No rows of data were given.")
  }

  if (grepl("chr", bed_df$chr)[1]){
    bed_df <- bed_df %>%
      dplyr::mutate(chr = paste0("chr", chr))
  }

  if (any(grepl("23", bed_df$chr))){
    bed_df <- bed_df %>%
      dplyr::mutate(chr = gsub("23", "X", chr))
  }

  if (any(grepl("24", bed_df$chr))){
    bed_df <- bed_df %>%
      dplyr::mutate(chr = gsub("24", "Y", chr))
  }

  # Make GRanges for CNV data
  granges <- GenomicRanges::GRanges(
    seqnames = bed_df$chr,
    ranges = IRanges::IRanges(
      start = bed_df$start,
      end = bed_df$end
    )
  )

  return(granges)
}

bed_overlap <- function(bed_file_1,
                        bed_file_2,
                        name_1,
                        name_2,
                        plot_name) {

  # Read in the BED files as GenomicRanges objects
  bed_granges_1 <- bed_to_granges(bed_file_1)
  bed_granges_2 <- bed_to_granges(bed_file_2)

  # Reduce these ranges 
  bed_granges_1 <- GenomicRanges::reduce(bed_granges_1)
  bed_granges_2 <- GenomicRanges::reduce(bed_granges_2)

  # Find intersection
  overlaps <- GenomicRanges::intersect(bed_granges_1, bed_granges_2)

  # Percent of TCGA Target region covered by overlap
  sum(overlaps@ranges@width) / sum(bed_granges_2@ranges@width)

  # Percent of PBTA Target region covered by overlap
  sum(bed_granges_1@ranges@width)

  # Make filename to save plot as
  plot.file <- file.path(plots_dir, paste0(plot_name, ".png"))

  # Make the Venn diagram
  grid::grid.newpage();
  venn.plot <- VennDiagram::draw.pairwise.venn(
    area1 = sum(bed_granges_1@ranges@width),
    area2 = sum(bed_granges_2@ranges@width),
    cross.area = sum(overlaps@ranges@width),
    category = c(name_1, name_2),
    fill = c("#F8766D", "#00BFC4"),
    cex = 2,
    cat.cex = 1.5,
    #cat.dist = c(-0.05, -0.05),
    #cat.pos = c(170, 65),
    ext.pos = 0,
    ext.dist = -0.01,
    ext.length = .8,
    ext.line.lwd = 2,
    ext.line.lty = "dashed",
    margin = 0.1);
  grid::grid.draw(venn.plot)

  png(plot.file);
  grid::grid.draw(venn.plot);
  dev.off()

  cat(" Ratio of", name_1, "overlapped:",
     sum(overlaps@ranges@width)/sum(bed_granges_1@ranges@width), "\n",
     "Ratio of", name_2, "overlapped:",
     sum(overlaps@ranges@width)/sum(bed_granges_2@ranges@width), "\n"
    )
}
