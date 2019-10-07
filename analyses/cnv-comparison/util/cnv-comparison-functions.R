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

read_in_cnv <- function(input_directory, file_path){
  # Given the file path of the CNV data, read in the data.
  #
  # Arg:
  #   input_directory: file path of input directory
  #   file_path: file path of input data
  #
  # Return:
  #   dataframe: the data frame containing the input data

  # Read in cnv data
  dataframe <-
    read.table(gzfile(file.path(input_directory, file_path)), header = TRUE)

  # Rename the columns
  dataframe <-
    dataframe[c("chrom", "loc.start", "loc.end", "ID", "num.mark", "seg.mean")]

  return(dataframe)

}

filter_segmean <- function(dataframe){
  # Given the data.frame containing CNV caller output, return a data.frame with
  # a column with the values 0 and 1 to represent loss and gain, respectively.
  #
  # Args:
  #   dataframe: data.frame containing CNV caller output
  #
  # Return:
  #   cnv_matrix: the data.frame given, returned with a column containing an
  #               aberration label for loss or gain

  # Rename columns for GenVisR function downstream
  colnames(dataframe) <-
    c("chromosome",
      "start",
      "end",
      "sample",
      "probes",
      "segmean")

  cnv_matrix <- dataframe %>%
    dplyr::mutate(aberration = "NA") %>%
    dplyr::mutate(aberration = dplyr::case_when(segmean < -0.5 ~ 0,
                                  segmean > 0.5 ~ 1)) %>%
    dplyr::filter(!is.na(aberration))

  return(cnv_matrix)

}

plot_cnFreq <-
  function(filtered_dataframe,
           plot_type,
           plot_title) {
    # Given the data.frame filtered for size of aberrations, plot the proportion
    # or frequency (denoted by the plot_type argument) of aberrations across
    # chromosomes with `GenVisR::cnFreq`
    #
    # Args:
    #   filtered_dataframe: data.frame filtered for cutoff size of aberrations
    #   plot_type: the type of plot (proportion or frequency)
    #   plot_title: a title string for the plot produced
    #
    # Return:
    #   aberration_plot: plot depicting the aberrations detected across
    #                    chromosomes

    aberration_plot <- GenVisR::cnFreq(
      filtered_dataframe,
      genome = "hg38",
      CN_low_cutoff = 0,
      CN_high_cutoff = .2,
      plot_title = toupper(plot_title),
      plotType = plot_type
    )

    return(aberration_plot)

  }

plot_boxplot <- function(dataframe, plot_title) {
  # Given the data.frame filtered for size of aberrations, plot the proportion
  # or frequency (denoted by the plot_type argument) of aberrations across
  # chromosomes with `geom_boxplot`
  #
  # Args:
  #   dataframe: data.frame filtered for cutoff size of aberrations
  #   plot_title: a title string for the plot produced
  #
  # Return:
  #   boxplot: boxplot depicting the log2 transformed sizes of aberrations
  #            detected across chromosomes

  # Create boxplot where the y-axis represents the log2 transformed segmean
  boxplot <- ggplot2::ggplot(dataframe,
                             ggplot2::aes(x = chromosome,
                                          y = (log2(segmean) + 1))) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::ylab("Log transformed segmean values") +
    ggplot2::ggtitle(toupper(plot_title))

  return(boxplot)

}

plot_histology_barplot <- function(dataframe, metadata, plot_title) {
  # Given the data.frame filtered by cutoff segmean value, plot the proportion
  # of aberrations across chromosomes with `geom_barplot`
  #
  # Args:
  #   dataframe: data.frame filtered for cutoff size of aberrations and
  #              joined with the metadata
  #   metadata: the relevant metadata data.frame
  #   plot_title: a title string for the plot produced
  #
  # Return:
  #   barplot: barplot depicting the proportion of aberrations detected across
  #            chromosomes, annoated by the broad_histology variable in the
  #            metadata

  # Create a data.frame with the filtered dataframe joined with the metadata
  meta_joined <- dataframe %>%
    dplyr::inner_join(metadata, by = c("sample" = "Kids_First_Biospecimen_ID"))

  # Create barplot where the y-axis represents the size of aberration
  barplot <- ggplot2::ggplot(meta_joined,
                             ggplot2::aes(x = chromosome,
                                          y = broad_histology,
                                          fill = broad_histology)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(toupper(plot_title))

  return(barplot)

}

plot_aberration_barplot <- function(dataframe, plot_title) {
  # Given the data.frame filtered by cutoff segmean value, plot the proportion
  # of aberrations across chromosomes with `geom_barplot`
  #
  # Args:
  #   dataframe: data.frame filtered by cutoff segmean value
  #   plot_title: a title string for the plot produced
  #
  # Return:
  #   barplot: barplot depicting the proportion of aberrations detected across
  #            chromosomes


  # Create barplot where the y-axis represents the size of aberration
  barplot <- ggplot2::ggplot(dataframe,
                             ggplot2::aes(x = chromosome,
                                          y = as.factor(aberration),
                                          fill = as.factor(aberration))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::ylab("Proportion of Aberrations") +
    ggplot2::labs(fill = "Loss(0)/Gain(1)") +
    ggplot2::ggtitle(toupper(plot_title))

  return(barplot)

}

plot_cowplot <- function(plot_a, plot_b, output_path, plot_name){
  # Given two plots, create a combined cowplot and save as a PDF in the
  # specified directory
  # Args:
  #   plot_a: the first plot, which will be positioned on top
  #   plot_b: the second plot, which will be positioned below the first
  #   output_path: the file.path to the output directory
  #   plot_name: the name the plot should be saved as

  # Save a combined cowplot plot of the cnvkit and controlfreec plots
  grid <- cowplot::plot_grid(plot_a, plot_b, ncol = 1)

  cowplot::save_plot(
    file.path(output_path, plot_name),
    grid,
    base_height = 12,
    base_width = 30
  )

}
