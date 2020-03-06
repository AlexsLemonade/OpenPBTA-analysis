# This script defines custom functions to be sourced in the
# notebooks of this module.
#
# Chante Bethell for CCDL 2020
#
# # #### USAGE
# This script is intended to be sourced in the
# 'analyses/compare-gistic/01-GISTIC-cohort-vs-histology-comparison.Rmd' and
# 'analyses/compare-gistic/02-GISTIC-tidy-data-prep.Rmd` notebooks as
# follows:
#
# source(file.path("util", "GISTIC-comparison-functions.R"))

genes_file = lgat_amp_genes_file
#### Implemented in both `01-GISTIC-cohort-vs-histology.Rmd` and `02-GISTIC-tidy-data-prep.Rmd`
format_gistic_genes <- function(genes_file,
                                include_peak_info = TRUE) {
  # Given the GISTIC result `amp_genes_conf_90.txt` or `del_genes_conf_90.txt`
  # files for the entire cohort or a specific histology, get the vector/data.frame
  # of genes for each amplification/deletion.
  #
  # Args:
  #  genes_file: file path to the `amp_genes.conf_90.txt` or
  #              `del_genes.conf_90.txt` file
  #
  # Return:
  #  genes_output: a vector/data.frame with all the genes included in the file
  
  genes_ragged_list <- data.table::fread(genes_file,
                                         data.table = FALSE)
  
  # Transpose data
  genes_transposed <- genes_ragged_list %>%
    t()
  
  # Make the header information in the first row the column names
  colnames(genes_transposed) <- genes_transposed[1,]
  
  genes_df <- genes_transposed %>%
    # Make into a data.frame object
    as.data.frame(stringsAsFactors = FALSE) %>%
    # Remove the row with header information
    dplyr::filter(cytoband != "cytoband")
  
  # Gather the gene data
  genes_output <- genes_df %>%
    tidyr::gather("wide_peak",
                  "gene",
                  -cytoband,
                  -`q value`,
                  -`residual q value`,
                  -`wide peak boundaries`) %>%
    dplyr::select(-wide_peak) %>%
    # Remove the blanks in the `gene` column that result from the gathering step
    dplyr::filter(gene != "")
  
  if (include_peak_info == FALSE) {
    # Return only the genes information
    genes_output <- genes_output %>%
      dplyr::select(gene) 
  }
  return(genes_output)
}

#### Implemented in `01-GISTIC-cohort-vs-histology-comparison.Rmd` ------------

# Code adapted from `analyses/cnv-chrom-plot/gistic_plot.Rmd`
plot_gistic_scores <- function(gistic_scores_file) {
  # Given the file path to a `scores.gistic` file, plot the gistic scores.
  #
  # Args:
  #   gistic_scores_file: file path to `scores.gistic` file
  #
  # Return:
  #    A ggbio plot of the given gistic scores

  # Read in and format gistic scores data
  gistic_scores <- data.table::fread(gistic_scores_file,
    data.table = FALSE
  ) %>%
    dplyr::rename("gscore" = "G-score") %>%
    # Recode 23 and 24 as X and Y.
    dplyr::mutate(
      Chromosome = as.character(Chromosome),
      Chromosome = dplyr::recode(Chromosome,
        "23" = "X",
        "24" = "Y"
      ),
      # Turn `Del` scores into negative `G-scores`
      # This is how GISTIC shows the scores.
      gscore = dplyr::case_when(
        Type == "Del" ~ -gscore,
        TRUE ~ gscore
      )
    )

  # Make GISTIC data into GRanges object
  gistic_ranges <- GenomicRanges::GRanges(
    seqnames = gistic_scores$Chromosome,
    ranges = IRanges::IRanges(
      start = gistic_scores$Start,
      end = gistic_scores$End
    ),
    score = gistic_scores$gscore,
    mcols = gistic_scores
  )

  # Plot the GISTIC scores
  gistic_plot <-
    ggbio::autoplot(
      gistic_ranges,
      ggplot2::aes(y = score, fill = mcols.Type),
      geom = "bar",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      size = 3,
      angle = 45,
      hjust = 1
    )) +
    colorblindr::scale_fill_OkabeIto(name = "Type") +
    ggplot2::ylab("G-scores") +
    ggplot2::scale_y_continuous(limits = c(-1, 1.2), breaks = seq(-1, 1.2, 0.2))

  # Return plot
  return(gistic_plot@ggplot)
}

plot_venn_diagram <-
  function(cohort_genes_vector,
             histology_genes_vector,
             histology_label) {
    # Given the GISTIC result `amp_genes_conf_90.txt` or `del_genes_conf_90.txt`
    # files for the entire cohort and a specific histology, plot a Venn Diagram
    # showing the counts of rows that overlap/do not overlap between the two
    # files.
    #
    # Args:
    #   cohort_genes_vector: a vector of all the genes in the amp/del file for
    #                           the entire cohort
    #   histology_genes_vector: a vector of all the genes in the amp/del file
    #                              for an individual histology
    #   histology_label: string indicating the individual histology, for the
    #                    purpose of labeling the Venn Diagram plots
    #
    # Returns:
    #  This function displays a Venn Diagram that represents the data that
    #  overlaps/does not overlap between the two given files

    # Define list for input in `venn` function
    input <-
      list(cohort_genes_vector, histology_genes_vector)

    names(input) <- c("cohort_genes", histology_label)

    # Display Venn Diagram
    gplots::venn(input)
  }

plot_genes_venn_diagram_wrapper <- function(cohort_genes_file,
                                            lgat_genes_file,
                                            hgat_genes_file,
                                            medulloblastoma_genes_file) {
  # Given the GISTIC result `amp_genes_conf_90.txt` or `del_genes_conf_90.txt`
  # files for the entire cohort and the three individual histologies, run the
  # `format_gistic_genes` and `plot_venn_diagram` functions to plot the overlaps
  # between the results for the entire cohort and the each of the individual
  # histologies.
  #
  # Args:
  #    cohort_genes_file: file path to the `amp_genes.conf_90.txt` or
  #                       `del_genes.conf_90.txt` file for the entire
  #                       cohort
  #    lgat_genes_file: file path to the `amp_genes.conf_90.txt` or
  #                     `del_genes.conf_90.txt` file for the LGAT
  #                     histology
  #    hgat_genes_file: file path to the `amp_genes.conf_90.txt` or
  #                     `del_genes.conf_90.txt` file for the HGAT
  #                     histology
  #    medulloblastoma_genes_file: file path to the `amp_genes.conf_90.txt` or
  #                                `del_genes.conf_90.txt` file for the
  #                                 medulloblastoma histology

  # Run `format_gistic_genes` function on each of the files
  cohort_genes_vector <- format_gistic_genes(cohort_genes_file, include_peak_info = FALSE)
  lgat_genes_vector <- format_gistic_genes(lgat_genes_file, include_peak_info = FALSE)
  hgat_genes_vector <- format_gistic_genes(hgat_genes_file, include_peak_info = FALSE)
  medulloblastoma_genes_vector <- format_gistic_genes(medulloblastoma_genes_file, include_peak_info = FALSE)

  # Run `plot_venn_diagram` for each comparison case
  lgat_venn <- plot_venn_diagram(cohort_genes_vector, lgat_genes_vector, "lgat_genes")
  hgat_venn <- plot_venn_diagram(cohort_genes_vector, hgat_genes_vector, "hgat_genes")
  medulloblastoma_venn <- plot_venn_diagram(cohort_genes_vector, medulloblastoma_genes_vector, "medulloblastoma_genes")

  # Save plots to list
  venn_plot_list <- list(lgat_venn, hgat_venn, medulloblastoma_venn)

  # Return the plot list
  return(venn_plot_list)
}

#### Implemented in `02-GISTIC-tidy-data-prep.Rmd` ----------------------------
all_lesions_file = cohort_all_lesions_file
amp_genes_file = cohort_amp_genes_file
del_genes_file = cohort_del_genes_file

prepare_gene_level_gistic <- function(all_lesions_file,
                                      amp_genes_file,
                                      del_genes_file,
                                      output_filename) {

  # Given the file paths to GISTIC's `all_lesion.conf_90.txt`,
  # `amp_genes.conf_90.txt`, and `del_genes.conf_90.txt` files,
  # read in and tidy this data into a data.frame that contains
  # the amplified/deleted genes, the detection peak, and
  # the IDs of the samples in the corresponding peaks along with
  # with their CN status call.
  #
  # Args:
  #   all_lesions_file: file path to GISTIC's `all_lesion.conf_90.txt` file
  #   amp_genes_file: file path to GISTIC's `amp_genes.conf_90.txt` file
  #   del_genes_file: file path to GISTIC's `del_genes.conf_90.txt` file
  #   output_filename: string to save the output file as
  #
  # Return:
  #   final_df: data.frame with the relevant data from the `all_lesions`,
  #             `amp_genes`, and `del_genes` files

  # Read in `all_lesions_file`
  gistic_all_lesions_df <- data.table::fread(all_lesions_file,
    data.table = FALSE
  ) %>% 
    # One half of the data is a copy of the other and the Unique Names
    # with `-CN values` have the actual copy number values. 
    dplyr::filter(grepl("- CN values", `Unique Name`))

  # Run `format_gistic_genes` function on `amp_genes` and `del_genes` files to get
  # a data.frame with just the genes and their corresponding detection peak
  amp_genes_df <- format_gistic_genes(amp_genes_file)

  del_genes_df <- format_gistic_genes(del_genes_file)
  
  # Bind the rows from the above data.frames into one data.frame
  final_df <- bind_rows("amp" = amp_genes_df, 
                        "del" = del_genes_df, 
                        .id = "direction")
  
  # Wrangle the GISTIC all lesions data to be in a comparable format with our CN calls
  gistic_all_lesions_df <- gistic_all_lesions_df %>%
    dplyr::select(dplyr::starts_with("BS"), 
                  "Unique Name", 
                  "Wide Peak Limits"
    ) %>%
    tidyr::gather(
      "Kids_First_Biospecimen_ID",
      "status",
      -`Unique Name`,
      -`Wide Peak Limits`
    ) %>%
    arrange(desc(`Wide Peak Limits`)) %>%
    dplyr::mutate(status = dplyr::case_when(
      status < 0 ~ "loss",
      status > 0 ~ "gain",
      status == 0 ~ "neutral"
    )) %>% 
    dplyr::mutate(
      # Keep just the chromosomal coordinates in the `Wide Peak Limits` column
      # (In order to merge with the amp/del genes data.frame)
      `Wide Peak Limits` = gsub("[(].*", "", `Wide Peak Limits`)
      )

  # Merge the data from the `all_lesions.conf_90.txt` file with the amp/del gene
  # data.frame prepped above
  final_df <- final_df %>%
    dplyr::left_join(gistic_all_lesions_df, by = c("wide peak boundaries" = "Wide Peak Limits")) %>%
    dplyr::select(gene, Kids_First_Biospecimen_ID, status, detection_peak = `Unique Name`)
  
  # Save data.frame to file
  readr::write_tsv(
    final_df,
    output_filepath
  )
}

prepare_cytoband_level_gistic <- function(all_lesions_file,
                                          output_filename) {

  # Given the file path to GISTIC's `all_lesion.conf_90.txt` file,
  # read in and tidy this data into a data.frame that contains
  # the cytoband, sample IDs, and their corresponding CN status call.
  #
  # Args:
  #   all_lesions_file: file path to GISTIC's `all_lesion.conf_90.txt` file
  #   output_filename: string to save the output file as
  #
  # Return:
  #   final_df: data.frame with the relevant data from the `all_lesions` file

  # Read in files
  gistic_all_lesions_df <- data.table::fread(all_lesions_file,
    data.table = FALSE
  )


  # Wrangle the GISTIC all lesions data to be in a comparable format with our CN calls
  gistic_all_lesions_df <- gistic_all_lesions_df %>%
    dplyr::select(
      -c(
        "Unique Name",
        "Wide Peak Limits",
        "Peak Limits",
        "Region Limits",
        "q values",
        "Residual q values after removing segments shared with higher peaks",
        "Broad or Focal",
        "Amplitude Threshold"
      )
    ) %>%
    as.data.frame() %>%
    tidyr::gather(
      "Kids_First_Biospecimen_ID",
      "status",
      -Descriptor
    ) %>%
    dplyr::mutate(status = dplyr::case_when(
      status < 0 ~ "loss",
      status > 0 ~ "gain",
      status == 0 ~ "neutral"
    )) %>%
    # The `select` function above got rid of some extra fields (fields that are
    # not needed for this analysis) from the `gistic_all_lesions_df`
    # object -- `distinct` removes any duplicate rows resulting from
    # the removal of the extra variables
    dplyr::distinct() %>%
    dplyr::rename(cytoband = Descriptor)

  # Save data.frame to file (for later cytoband level comparison)
  readr::write_tsv(
    gistic_all_lesions_df,
    file.path(results_dir, output_filename)
  )
}
