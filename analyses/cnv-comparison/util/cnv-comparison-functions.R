# This script defines filtering and plotting functions to be sourced in 
# `01-cnv-comparison-plotting.R`
#
# Chante Bethell for CCDL 2019
# 
# Usage:
# This script is intended to be run via the command line from the top directory 
# of the repository as follows:
#
# Rscript analyses/cnv-comparison/util/cnv-comparison-functions.R

read_in_cnv <- function(file_path){
  # Given the file path of the CNV data, read in the data.
  #
  # Arg:
  #   file_path: file path of input data
  #
  # Return:
  #   cnv_df: the data frame containing the input data

  # Read in cnv data
  cnv_df <-
    read.table(gzfile(file_path), header = TRUE, stringsAsFactors = FALSE)

  # Rearrange the columns
  cnv_df <- cnv_df %>%
    dplyr::select("chrom", "loc.start", "loc.end", "ID", "num.mark", "seg.mean") %>%
    dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y"))))

  return(cnv_df)

}

filter_segmean <- function(cnv_df, segmean_cutoff){
  # Given the data.frame containing CNV caller output, return a data.frame with
  # a column with the values 0 and 1 to represent loss and gain, respectively.
  #
  # Args:
  #   cnv_df: data.frame containing CNV caller output
  #   segmean_cutoff: segment mean value cutoff 
  #
  # Return:
  #   filtered_cnv_df: the data.frame given, returned with a column containing
  #                    an aberration label for loss or gain

  # Rename columns for GenVisR function downstream
  colnames(cnv_df) <-
    c("chromosome",
      "start",
      "end",
      "sample",
      "probes",
      "segmean")

  # Define an aberration column 
  filtered_cnv_df <- cnv_df %>%
    dplyr::mutate(aberration = dplyr::case_when(segmean < -segmean_cutoff ~ 0,
                                  segmean > segmean_cutoff ~ 1)) %>%
    dplyr::filter(!is.na(aberration))

  filtered_cnv_df$aberration <- 
    factor(filtered_cnv_df$aberration, labels = c("Loss", "Gain"))
  
  return(filtered_cnv_df)

}

plot_violin <- function(filtered_cnv_df) {
  # Given the data.frame filtered for size of aberrations, plot the proportion
  # or frequency (denoted by the plot_type argument) of aberrations across
  # chromosomes with `geom_violin`
  #
  # Args:
  #   filtered_cnv_df: data.frame filtered for cutoff size of aberrations
  #
  # Return:
  #   violin_plot: violin plot depicting the log2 transformed sizes of
  #                aberrations detected across chromosomes

  # Create violin plot where the y-axis represents the log2 transformed segmean
  violin_plot <- ggplot2::ggplot(filtered_cnv_df,
                             ggplot2::aes(x = chromosome,
                                          y = (log2(abs(segmean)) + 1))) +
    ggplot2::geom_violin() +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(cnv_caller ~ .) +
    ggplot2::ylab("Log transformed segmean values") +
    ggplot2::theme(strip.text.y = ggplot2::element_text(size = 14))
    
  return(violin_plot)

}

plot_histology_barplot <- function(filtered_cnv_df, metadata) {
  # Given the data.frame filtered by cutoff segmean value with `filter_segmean`,
  # plot the proportion of aberrations across chromosomes with `geom_barplot`
  #
  # Args:
  #   filtered_cnv_df: data.frame filtered for cutoff size of aberrations and
  #                    joined with the metadata
  #   metadata: the relevant metadata data.frame
  #
  # Return:
  #   barplot: barplot depicting the proportion of aberrations detected across
  #            chromosomes, annotated by the `broad_histology`` variable in the
  #            metadata

  # Create a data.frame with the filtered dataframe joined with the metadata
  meta_joined <- filtered_cnv_df %>%
    dplyr::inner_join(metadata, by = c("sample" = "Kids_First_Biospecimen_ID"))

  # Create barplot where the y-axis represents the size of aberration
  barplot <- ggplot2::ggplot(meta_joined,
                             ggplot2::aes(x = chromosome,
                                          fill = broad_histology)) +
    ggplot2::geom_bar() +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(cnv_caller ~ aberration) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y = ggplot2::element_text(size = 14))
  
  return(barplot)

}

plot_aberration_barplot <- function(filtered_cnv_df) {
  # Given the data.frame filtered by cutoff segmean value with `filter_segmean`,
  # plot the proportion of aberrations across chromosomes with `geom_barplot`
  #
  # Args:
  #   filtered_cnv_df: data.frame filtered by cutoff segmean value
  #
  # Return:
  #   barplot: barplot depicting the proportion of aberrations detected across
  #            chromosomes
  
  # Create barplot where the y-axis represents the size of aberration
  barplot <- ggplot2::ggplot(filtered_cnv_df,
                             ggplot2::aes(x = chromosome,
                                          fill = aberration)) +
    ggplot2::geom_bar() +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(cnv_caller ~ aberration) +
    ggplot2::ylab("Proportion of Aberrations") +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y = ggplot2::element_text(size = 14))

  return(barplot)

}

plot_cowplot <- function(plot_list, output_path, plot_name){
  # Given two plots, create a combined cowplot and save as a PDF in the
  # specified directory
  # Args:
  #   plot_list: list of plots
  #   output_path: the file.path to the output directory
  #   plot_name: the name the plot should be saved as

  # Save a combined cowplot plot of the cnvkit and controlfreec plots
  grid <-
    cowplot::plot_grid(
      plot_list[[1]],
      plot_list[[2]],
      ncol = 1,
      labels = c("CNVkit", "ControlFREEC"),
      label_y = 1.05
    ) + ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 2), "cm"))

  cowplot::save_plot(
    file.path(output_path, plot_name),
    grid,
    base_height = 12,
    base_width = 30
  )
  
  return(grid)
  
}
