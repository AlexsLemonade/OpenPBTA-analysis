# Intended for import only
# Chante Bethell for CCDL 2019
#
# Assign function to perform the dimension reduction techniques
perform_dimension_reduction <- function(transposed_expression_matrix,
                                        method,
                                        model_filename,
                                        output_directory,
                                        perplexity_parameter,
                                        neighbors_parameter,
                                        seed = 2019) {
  # Given a data.frame that containes the transposed expression matrix and the
  # name of a dimension reduction method, perform the dimension reduction
  # technique on the transposed matrix
  #
  # Args:
  #   transposed_expression_matrix: Transposed `RSEM` or `Kallisto` expression
  #                                 matrix
  #   method: Dimension reduction method to be performed
  #   model_filename: Name of the output model file
  #   output_directory: file.path to the output directory
  #   perplexity_parameter: integer defining the perplexity parameter for t-SNE
  #   neighbors_parameter: integer defining the n_neighbors parameter for UMAP
  #   seed: Seed, set to a default of 2019
  #
  # Returns:
  #   dimension_reduction_df: data.frame containing the resulting scores of the
  #                           dimension reduction technique, with a column that
  #                           includes sample identifiers

  # Set seed for reproducibility
  set.seed(seed)

  # Save rownames as a vector
  ID <- rownames(transposed_expression_matrix)

  # Perform dimension reduction
  if (method == "PCA") {
    dimension_reduction_results <- prcomp(transposed_expression_matrix)
    dimension_reduction_df <-
      data.frame(dimension_reduction_results$x)
  } else if (method == "t-SNE") {
    dimension_reduction_results <-
      Rtsne::Rtsne(transposed_expression_matrix,
                   perplexity = perplexity_parameter
      )
    dimension_reduction_df <-
      data.frame(dimension_reduction_results$Y)
  } else if (method == "UMAP") {
    dimension_reduction_results <-
      umap::umap(transposed_expression_matrix,
                 n_neighbors = neighbors_parameter
      )
    dimension_reduction_df <-
      data.frame(dimension_reduction_results$layout)
  } else {
    stop(paste(method, "is not a supported argument"))
  }

  # Assign the relevant rownames(Kids_First_Biospecimen_ID) to a column
  dimension_reduction_df$Kids_First_Biospecimen_ID <- ID

  # Write the resulting model to a rds file
  readr::write_rds(
    dimension_reduction_results,
    file.path(
      output_directory,
      paste(model_filename, "_model.rds", sep = "")
    )
  )

  return(dimension_reduction_df)
}

# Assign function to align the metadata with the dimension reduction scores
align_metadata <- function(dimension_reduction_df,
                           metadata_df,
                           scores_filename,
                           output_directory) {
  # Given a data.frame that contains scores from a dimensionality reduction
  # technique, align the metadata to the data.frame in preparation for plotting
  #
  # Note: The rownames provided in the id argument are synonymous with the
  # metadata's `Kids_First_Biospecimen_ID`
  #
  # Args:
  #   dimension_reduction_df: Name of the data.frame containing dimension
  #   reduction scores
  #   metadata_df : Name of the data.frame containing the relevant metadata
  #   scores_filename: Name of the output scores tsv file
  #   output_directory: file.path to the output directory
  #   strategy: The selection strategy used for RNA Sequencing
  #
  # Returns:
  #   aligned_scores_df: A data.frame containing the dimension reduction scores
  #         along with the variables from `metadata_df`

  # Join the dimension reductions scores data.frame and the metadata data.frame
  aligned_scores_df <- dimension_reduction_df %>%
    dplyr::inner_join(metadata_df, by = c("Kids_First_Biospecimen_ID"))

  # Write the resulting metadata aligned data.frame to a tsv file
  readr::write_tsv(
    aligned_scores_df,
    file.path(
      output_directory,
      paste(scores_filename, "_scores_aligned.tsv", sep = "")
    )
  )

  return(aligned_scores_df)
}

# Assign wrapper function to execute the two functions above
dimension_reduction_wrapper <- function(transposed_expression_matrix,
                                        method,
                                        seed,
                                        metadata_df,
                                        filename_lead,
                                        output_directory,
                                        perplexity_parameter = NULL,
                                        neighbors_parameter = NULL) {
  # Given a transposed expression matrix, the metadata data.frame, the name of
  # the dimension reduction method to be performed and a seed with a default
  # value, perform the dimension reduction, and align the metadata with the
  # resulting scores
  #
  # Args:
  #   transposed_expression_matrix: Transposed `RSEM` or `Kallisto` expression
  #                                 matrix
  #   method: Dimension reduction method to be performed
  #   seed: seed set to default 2019
  #   metadata_df: Name of the data.frame containing the relevant metadata
  #   filename_lead: Character string that will be used to name the output files
  #   output_directory: file.path to the output directory
  #   perplexity_parameter: integer defining the perplexity parameter for t-SNE
  #   neighbors_parameter: integer defining the n_neighbors parameter for UMAP
  #
  # Returns:
  #   aligned_scores_df: data.frame containing dimension reduction scores and
  #                      aligned metadata variables
  #

  # Run `perform_dimension_reduction` function
  dimension_reduction_df <-
    perform_dimension_reduction(
      transposed_expression_matrix,
      method,
      model_filename = filename_lead,
      output_directory,
      perplexity_parameter,
      neighbors_parameter,
      seed = seed
    )

  # Run `align_metadata` function
  aligned_scores_df <-
    align_metadata(
      dimension_reduction_df,
      metadata_df = metadata_df,
      scores_filename = filename_lead,
      output_directory
    )

}

# Function to plot
plot_dimension_reduction <- function(aligned_scores_df,
                                     point_color,
                                     point_shape = NULL,
                                     point_size = NULL,
                                     x_label,
                                     y_label,
                                     alpha_value = 0.3,
                                     score1 = 1,
                                     score2 = 2,
                                     color_palette = NULL) {
    # Given a data.frame that contains the scores of a dimension reduction
    # analysis and the information we want from the metadata, make a scatterplot.
    #
    # Args:
    #   aligned_scores_df: data.frame containing dimension reduction scores,
    #                           and relative information from the metadata
    #   point_color: the variable whose information will be used to color
    #                the points on the plot
    #   point_shape: the variable whose information will be used to shape
    #                the points on the plot, this is NULL by default
    #   point_size: the variable whose information will be used to size the
    #               points on the plot, this is NULL by default
    #   x_label: the x-axis label, character
    #   y_label: the y-axis label, character
    #   score1: the column number of the first dimension reduction score that we
    #           want to plot, 1 by default
    #   score2: the column number of the second dimension reduction score that
    #           we want to plot, 2 by default
    #   color_palette: color palette to be used for the point color - one column
    #                  should contain the same values as what is in point_color
    #                  and the second column should contain hex codes
    #
    # Returns:
    #   dimension_reduction_plot: the plot representing the dimension reduction
    #                             scores in the given data.frame, using the
    #                             values in `point_color` as symbols to color
    #                             the points

    if (!(point_color %in% colnames(aligned_scores_df))) {
      stop(paste(point_color, "is not column in aligned_scores_df"))
    }

    if (!(is.null(point_shape)) && !(point_shape %in% colnames(aligned_scores_df))) {
      stop(paste(point_shape, "is not column in aligned_scores_df"))
    }

    if (!(is.null(point_size)) && !(point_size %in% colnames(aligned_scores_df))) {
      stop(paste(point_size, "is not column in aligned_scores_df"))
    }

    # transform the strings `point_color` into symbols for plotting
    color_sym <- rlang::sym(point_color)

    # transform the strings `point_size` and `point_shape` into symbols for
    # plotting when they are not NULL
    if (!(is.null(point_size))) {
      point_size <- rlang::sym(point_size)
    }

    if (!is.null(point_shape)) {
      point_shape <- rlang::sym(point_shape)
    }

    `%>%` <- dplyr::`%>%`

    if (!is.null(color_palette)) {

      dimension_reduction_plot <- ggplot2::ggplot(
        aligned_scores_df,
        ggplot2::aes(
          x = dplyr::pull(aligned_scores_df, score1),
          y = dplyr::pull(aligned_scores_df, score2),
          color = !!color_sym,
          size = !!point_size,
          shape = !!point_shape
        )
      ) +
      ggplot2::scale_color_manual(values = color_palette)
    } else {
      dimension_reduction_plot <- ggplot2::ggplot(
        aligned_scores_df,
        ggplot2::aes(
          x = dplyr::pull(aligned_scores_df, score1),
          y = dplyr::pull(aligned_scores_df, score2),
          color = !!color_sym,
          size = !!point_size,
          shape = !!point_shape
        )
      ) +
      ggplot2::scale_color_manual(values = (c(
        "#2b3fff", "#ec102f", "#235e31",
        "#1a6587", "#11e38c", "#a22f80",
        "#fe5900", "#1945c5", "#51f310",
        "#8b20d3", "#799d10", "#881c23",
        "#3fc6f8", "#fe5cde", "#0a7fb2",
        "#f2945a", "#6b4472", "#f4d403",
        "#76480d", "#a6b6f9", "#0d76ff",
        "#f3e500", "#a72a1b", "#319106",
        "#0051bc", "#000000"
      )))
    }

    dimension_reduction_plot <- dimension_reduction_plot +
      ggplot2::geom_point(alpha = alpha_value) +
      ggplot2::labs(x = x_label, y = y_label) +
      ggplot2::theme_bw()

    return(dimension_reduction_plot)
}
