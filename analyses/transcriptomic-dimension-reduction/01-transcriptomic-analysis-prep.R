# Unsupervised Analysis of Transcriptomic Differences - Data Prep
# Chante Bethell for CCDL 2019
#
# This notebook is the first part of an analysis that demonstrates samples in the
# PBTA cluster by cancer type using dimensionality reduction techniques, namely,
# Principal Component Analysis (PCA),
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
# This script is intended to be run via the command line from the top directory of the repository as follows:
#
# Rscript --vanilla analyses/transcriptomic-dimension-reduction/01-transcriptomic-analysis-prep.R --perplexity 10 --neighbors 15
#
# where the integers can be replaced with other integers to change the 
# perplexity paramater for t-SNE and the neighbors parameter for UMAP. 


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
# magrittr pipe
`%>%` <- dplyr::`%>%`

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

# Assign function to perform the dimension reduction techniques
perform_dimension_reduction <-
  function(transposed_expression_matrix,
           method,
           seed,
           filename,
           output_directory,
           perplexity_parameter,
           neighbors_parameter) {
    # Given a data.frame that containes the transposed expression matrix and the
    # name of a dimension reduction method, perform the dimension reduction
    # technique on the transposed matrix
    # Args:
    #   transposed_expression_matrix: Transposed `RSEM` or `Kallisto` expression
    #                                 matrix
    #   method: Dimension reduction method to be performed
    #   seed: Seed set to a default of 2019
    #   filename: Name of the output file
    #   output_directory: file.path to the output directory
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
      prcomp_results <- prcomp(transposed_expression_matrix)
      dimension_reduction_df <- data.frame(prcomp_results$x[, 1:2])
    } else if (method == "t-SNE") {
      tsne_results <-
        Rtsne::Rtsne(transposed_expression_matrix, perplexity = perplexity_parameter)
      dimension_reduction_df <- data.frame(tsne_results$Y)
    } else if (method == "UMAP") {
      umap_results <- umap::umap(transposed_expression_matrix, n_neighbors = neighbors_parameter)
      dimension_reduction_df <- data.frame(umap_results$layout)
    } else {
      stop("method is not a supported argument")
    }

    # Assign the relevant rownames(Kids_First_Biospecimen_ID) to a column
    dimension_reduction_df$Kids_First_Biospecimen_ID <- ID
    # Write the resulting model to a rds file
    readr::write_rds(
      dimension_reduction_df,
      file.path(
        output_directory,
        "models",
        paste(filename, "_model.rds", sep = "")
      )
    )

    return(dimension_reduction_df)
    
  }
# Assign function to align the metadata with the dimension reduction scores
align_metadata <-
  function(dimension_reduction_df,
             metadata_df,
             filename,
             output_directory) {
    # Given a data.frame that contains scores from a dimensionality reduction
    # technique, align the metadata to the data.frame in preparation for plotting
    #
    # Note: The rownames provided in the id argument are synonymous with the
    # metadata's `Kids_First_Biospecimen_ID`
    #
    # Args:
    #   dimension_reduction_df: Name of the data.frame containing dimension reduction scores
    #   metadata_df : Name of the data.frame containing the relevant metadata
    #   filename: Name of the output file
    #   output_directory: file.path to the output directory
    #
    # Returns:
    #   aligned_scores_df: A data.frame containing the dimension reduction scores along with
    #         the variables from `metadata_df`

    # Join the dimension reductions scores data.frame and the metadata data.frame
    aligned_scores_df <- dimension_reduction_df %>%
      dplyr::inner_join(metadata_df, by = c("Kids_First_Biospecimen_ID"))

    # Write the resulting metadata aligned data.frame to a tsv file
    readr::write_tsv(
      aligned_scores_df,
      file.path(
        output_directory,
        "results",
        paste(filename, "_scores_aligned.tsv", sep = "")
      )
    )

    return(aligned_scores_df)
    
  }
# Assign wrapper function to execute the two functions above
dimension_reduction_wrapper <-
  function(transposed_expression_matrix,
           method,
           seed,
           metadata_df,
           filename,
           output_directory,
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
    #   filename: Filename for the RDS output file
    #   output_directory: file.path to the output directory
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
        seed = seed,
        filename,
        output_directory,
        perplexity_parameter,
        neighbors_parameter
      )
    # Run `align_metadata` function
    aligned_scores_df <-
      align_metadata(
        dimension_reduction_df,
        metadata_df,
        filename,
        output_directory
      )

    return(aligned_scores_df)
    
  }

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
output_dir <- file.path(root_dir, "analyses", "transcriptomic-dimension-reduction")
models_dir <- file.path(output_dir, "models")
results_dir <- file.path(output_dir, "results")

# Create directories to hold the output.
if (!dir.exists(models_dir)) {
  dir.create(models_dir)
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in metadata
metadata_df <-
  data.frame(readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")))

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

# RSEM
# Run `dimension_reduction_wrapper` function using PCA
rsem_pca_aligned_scores <-
  dimension_reduction_wrapper(transposed_rsem_data,
    "PCA",
    seed,
    metadata_df,
    "rsem_pca",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using t-SNE
rsem_tsne_aligned_scores <-
  dimension_reduction_wrapper(transposed_rsem_data,
    "t-SNE",
    seed,
    metadata_df,
    "rsem_tsne",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using UMAP
rsem_umap_aligned_scores <-
  dimension_reduction_wrapper(transposed_rsem_data,
    "UMAP",
    seed,
    metadata_df,
    "rsem_umap",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )

# Kallisto
# Run `dimension_reduction_wrapper` function using PCA
kallisto_pca_aligned_scores <-
  dimension_reduction_wrapper(transposed_kallisto_data,
    "PCA",
    seed,
    metadata_df,
    "kallisto_pca",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using t-SNE
kallisto_tsne_aligned_scores <-
  dimension_reduction_wrapper(transposed_kallisto_data,
    "t-SNE",
    seed,
    metadata_df,
    "kallisto_tsne",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )
# Run `dimension_reduction_wrapper` function using UMAP
kallisto_umap_aligned_scores <-
  dimension_reduction_wrapper(transposed_kallisto_data,
    "UMAP",
    seed,
    metadata_df,
    "kallisto_umap",
    output_dir,
    perplexity_parameter,
    neighbors_parameter
  )
