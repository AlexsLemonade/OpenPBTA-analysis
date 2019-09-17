# Unsupervised Analysis of Transcriptomic Differences - Data Prep
# Chante Bethell for CCDL 2019
#
# This script is the first part of an analysis that demonstrates samples in
# the PBTA cluster by cancer type using dimensionality reduction techniques,
# namely, Principal Component Analysis (PCA),
# t-Distributed Stochastic Neighbor Embedding (t-SNE),
# and Uniform Manifold Approximation and Projection (UMAP).
#
# This script focuses on the preparation of the expression data and the
# execution of the dimension reduction techniques.
#
# It addresses issue #9 (https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/9)
# in the Open-PBTA analysis repository.
#
# The output files include dimension reduction score models of the RSEM and
# Kallisto data for each of the three dimension reduction techniques mentioned
# above, as well as tsv files containing these scores aligned with the metadata.
#
# # Usage
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript --vanilla --max-ppsize=500000 analyses/transcriptomic-dimension-reduction/01-transcriptomic-analysis-prep.R --perplexity 10 --neighbors 15
#
# where the integers can be replaced with other integers to change the
# perplexity paramater for t-SNE and the neighbors parameter for UMAP.

#### Install packages ----------------------------------------------------------

# This is needed for running the t-SNE analysis
if (!("Rtsne" %in% installed.packages())) {
  install.packages("Rtsne")
}

# This is needed for running the umap analysis
if (!("umap" %in% installed.packages())) {
  install.packages("umap")
}

# This is needed for taking arguments from the command line
if (!("optparse" %in% installed.packages())) {
  install.packages("optparse")
}

#### Functions -----------------------------------------------------------------

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Assign function to perform the dimension reduction techniques
perform_dimension_reduction <- function(transposed_expression_matrix,
                                        method,
                                        seed,
                                        model_filename,
                                        output_directory,
                                        perplexity_parameter,
                                        neighbors_parameter,
                                        metadata_df,
                                        strategy) {
  # Given a data.frame that containes the transposed expression matrix and the
  # name of a dimension reduction method, perform the dimension reduction
  # technique on the transposed matrix
  # Args:
  #   transposed_expression_matrix: Transposed `RSEM` or `Kallisto` expression
  #                                 matrix
  #   method: Dimension reduction method to be performed
  #   seed: Seed set to a default of 2019
  #   model_filename: Name of the output model file
  #   output_directory: file.path to the output directory
  #   perplexity_parameter: integer defining the perplexity parameter for t-SNE
  #   neighbors_parameter: integer defining the n_neighbors parameter for UMAP
  #
  # Returns:
  #   dimension_reduction_df: data.frame containing the resulting scores of the
  #                           dimension reduction technique, with a column that
  #                           includes sample identifiers

  # Set seed for reproducibility
  set.seed(seed)
  
  # Filter for selection strategy 
  metadata_df2 <- metadata_df %>%
    dplyr::filter(RNA_library %in% strategy)
  
  transposed_expression_matrix <- transposed_expression_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% metadata_df2$Kids_First_Biospecimen_ID)
  
  transposed_expression_matrix <- transposed_expression_matrix %>%
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
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
    stop("method is not a supported argument")
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
  
  dimension_list <- (list(dimension_reduction_df, metadata_df2))
  
  return(dimension_list)
}

# Assign function to align the metadata with the dimension reduction scores
align_metadata <- function(dimension_list,
                           scores_filename,
                           output_directory,
                           strategy) {
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
  aligned_scores_df <- dimension_list[[1]] %>%
    dplyr::inner_join(dimension_list[[2]], by = c("Kids_First_Biospecimen_ID")) 

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
                                        model_filename,
                                        scores_filename,
                                        output_directory,
                                        strategy,
                                        perplexity_parameter,
                                        neighbors_parameter) {
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
  #   model_filename: Filename for the model RDS output file
  #   scores_filename: Filename for the scores tsv output file
  #   output_directory: file.path to the output directory
  #   strategy: The selection strategy used for RNA Sequencing
  #   perplexity_parameter: integer defining the perplexity parameter for t-SNE
  #   neighbors_parameter: integer defining the n_neighbors parameter for UMAP
  #
  # Returns:
  #   aligned_scores_df: data.frame containing dimension reduction scores and
  #                      aligned metadata variables
  #

  # Run `perform_dimension_reduction` function
  dimension_list <-
    perform_dimension_reduction(
      transposed_expression_matrix,
      method,
      seed,
      model_filename,
      output_directory,
      perplexity_parameter,
      neighbors_parameter,
      metadata_df,
      strategy
    )

  # Run `align_metadata` function
  aligned_scores_df <-
    align_metadata(
      dimension_list,
      scores_filename,
      output_directory,
      strategy
    )

  return(aligned_scores_df)
}

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-p", "--perplexity"),
    type = "integer",
    default = 10,
    help = "t-SNE perplexity parameter integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2019,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-n", "--neighbors"),
    type = "integer",
    default = 15,
    help = "UMAP neighbors parameter integer",
    metavar = "integer"
  )
)
# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Define perplexity parameter as the parameter passed via command line
perplexity_parameter <- opt$perplexity

# Define seed as the parmeter passed via the command line
seed <- opt$seed

# Define neighbors parameter as the parameter passed via command line
neighbors_parameter <- opt$neighbors

#### Directories ---------------------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directory
output_dir <- file.path(
  root_dir, "analyses",
  "transcriptomic-dimension-reduction", "results"
)

# Create directory to hold the output.
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata_df <-
  data.frame(readr::read_tsv(file.path(
    root_dir, "data",
    "pbta-histologies.tsv"
  )))

# Read in kallisto expression data
exp_kallisto <- data.frame(readr::read_rds(file.path(
  root_dir, "data", "pbta-gene-expression-kallisto.rds"
)))

# Read in RSEM expression data
exp_rsem <- data.frame(readr::read_rds(file.path(
  root_dir, "data", "pbta-gene-expression-rsem.fpkm.rds"
)))

# Transform the non-numeric "gene_id" column into rownames
exp_kallisto <- exp_kallisto[, -1] %>%
  dplyr::filter(!duplicated(gene_id)) %>%
  dplyr::filter(complete.cases(.)) %>%
  tibble::column_to_rownames("gene_id")

exp_rsem <- exp_rsem %>%
  dplyr::filter(complete.cases(.)) %>%
  tibble::column_to_rownames("gene_id")

# Filter out low count genes
gene_count_threshold <- 100

kallisto_genes_to_keep <-
  rowSums(exp_kallisto) >= gene_count_threshold
exp_kallisto <- exp_kallisto[kallisto_genes_to_keep, ]

rsem_genes_to_keep <- rowSums(exp_rsem) >= gene_count_threshold
exp_rsem <- exp_rsem[rsem_genes_to_keep, ]

# Transpose the data
transposed_rsem_data <- t(exp_rsem)
transposed_kallisto_data <- t(exp_kallisto)

#### RSEM ----------------------------------------------------------------------

# Run `dimension_reduction_wrapper` function using PCA
rsem_pca_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "PCA",
    seed,
    metadata_df,
    "rsem_pca",
    "polyA_rsem_pca",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_pca_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "PCA",
    seed,
    metadata_df,
    "rsem_pca",
    "stranded_rsem_pca",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_pca_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "PCA",
    seed,
    metadata_df,
    "rsem_pca",
    "rsem_pca",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )

# Run `dimension_reduction_wrapper` function using t-SNE
rsem_tsne_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "t-SNE",
    seed,
    metadata_df,
    "rsem_tsne",
    "polyA_rsem_tsne",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_tsne_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "t-SNE",
    seed,
    metadata_df,
    "rsem_tsne",
    "stranded_rsem_tsne",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_tsne_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "t-SNE",
    seed,
    metadata_df,
    "rsem_tsne",
    "rsem_tsne",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using UMAP
rsem_umap_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "UMAP",
    seed,
    metadata_df,
    "rsem_umap",
    "polyA_rsem_umap",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_umap_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "UMAP",
    seed,
    metadata_df,
    "rsem_umap",
    "stranded_rsem_umap",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
rsem_umap_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_rsem_data,
    "UMAP",
    seed,
    metadata_df,
    "rsem_umap",
    "rsem_umap",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )

#### Kallisto ------------------------------------------------------------------

# Run `dimension_reduction_wrapper` function using PCA
kallisto_pca_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "PCA",
    seed,
    metadata_df,
    "kallisto_pca",
    "polyA_kallisto_pca",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_pca_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "PCA",
    seed,
    metadata_df,
    "kallisto_pca",
    "stranded_kallisto_pca",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_pca_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "PCA",
    seed,
    metadata_df,
    "kallisto_pca",
    "kallisto_pca",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using t-SNE
kallisto_tsne_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "t-SNE",
    seed,
    metadata_df,
    "kallisto_tsne",
    "polyA_kallisto_tsne",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_tsne_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "t-SNE",
    seed,
    metadata_df,
    "kallisto_tsne",
    "stranded_kallisto_tsne",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_tsne_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "t-SNE",
    seed,
    metadata_df,
    "kallisto_tsne",
    "kallisto_tsne",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using UMAP
kallisto_umap_aligned_scores_polyA <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "UMAP",
    seed,
    metadata_df,
    "kallisto_umap",
    "polyA_kallisto_umap",
    output_dir,
    c("poly-A"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_umap_aligned_scores_stranded <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "UMAP",
    seed,
    metadata_df,
    "kallisto_umap",
    "stranded_kallisto_umap",
    output_dir,
    c("stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
kallisto_umap_aligned_scores <-
  dimension_reduction_wrapper(
    transposed_kallisto_data,
    "UMAP",
    seed,
    metadata_df,
    "kallisto_umap",
    "kallisto_umap",
    output_dir,
    c("poly-A", "stranded"),
    perplexity_parameter,
    neighbors_parameter
  )
