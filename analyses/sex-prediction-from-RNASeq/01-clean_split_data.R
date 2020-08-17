

#Functionality/purpose:

#This will clean the gene expression data -- i.e., it drops anything with invalid labels
#in any of the target columns specified and if the training partition size < 1.0,
#it will split the data into training and testing sets that can be saved as separate files to be used downstream.

#Inputs and arguments:
#gene expression matrix
#histologies file
#seed with some default
#training partition size with some default
#output filenames or output directory, depending on how it is implemented


#--------Library calls and run parameters.  Modify the assignments as needed for subsequent runs.

cat(paste("\n\nStart script 01_clean_split_data at", Sys.time(), "\n", sep=" "))

library(dplyr)
library(tibble)
library(readr)
library(optparse)

# We tune our predictive parameters using cross-validation on a training set and test the predictive accuracy
# of the optimal model on a holdout test set.
# train_percent represents the train/test partition; the seed is used to select training samples randomly.
# The test set is whatever did not get selected for the training set.
# train/test partition percentage and the seed are specified as command line arguments --train_percent and --seed and
# assigned to variables train_percent and seed below.

# Several gene expression files can serve as input to this analysis.
# The file name for this run is specified as command line argument --expression and assigned to gene_expression_file_name below.

# The output directory is specified as command line arument --output_directory and assigned to output_directory below.
# training and test data files are saved into the output directory.
# A TSV file that contains the biospecimen ID, and target columns
# information that has been cleaned is also saved to the output directory.

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--expression"),
    type = "character",
    default = NULL,
    help = "full path to expression data RDS file",
  ),
  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to metadata TSV file"
  ),
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory"
  ),
  optparse::make_option(
    c("-f", "--train_expression_file_name"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-a", "--test_expression_file_name"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-b", "--train_targets_file_name"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-c", "--test_targets_file_name"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-d", "--full_targets_file_name"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 36354,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-t", "--train_percent"),
    type = "double",
    default = 0.7,
    help = "The proportion of samples that should be in the training set when splitting into training and testing sets",
    metavar = "double"
  ),
  optparse::make_option(
    c("-e", "--target_columns"),
    type = "character",
    default = "reported_gender",
    help = "Names of columns to use as targets in testing separated with a comma",
  ),
  optparse::make_option(
    c("-r", "--train_target_column"),
    type = "character",
    default = "reported_gender",
    help = "Names of column in the target dataframe used as training label",
  )  
)


# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#set the run parameter variables from the command line arguments

train_percent <- opt$train_percent
seed <- opt$seed
set.seed(seed)

gene_expression_file_name <- opt$expression
histologies_file_name <- opt$metadata

target_column_list <- strsplit(opt$target_columns, ",")

# check train_target_column is in target_column_list
if (!(opt$train_target_column %in% target_column_list[[1]])) {
  stop(paste("TRAIN_TARGET_COLUMN value is not in target_Columns list. Check both arguments."))
}

# Create specified output directory if it does not yet exist
output_directory <- opt$output_directory
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}



#--------

#--------Read input data files

#RNA-Seq data is in pbta-gene-expression-kallisto.rds.

ge <- readRDS(gene_expression_file_name)


#Metadata is in pbta-histologies.tsv.  train_target_column is our target variable.
#Kids_First_Biospecimen_ID is the unique identifier for each sample.

histologies <- read.delim(histologies_file_name, header=TRUE, sep="\t", stringsAsFactors = FALSE)

samples_list <- list()

for (t in target_column_list[[1]]) {
  
  samples_list[[t]] <- histologies[which(toupper(histologies[, t]) == "MALE" |
                                           toupper(histologies[, t]) == "FEMALE"), "Kids_First_Biospecimen_ID"]
  
}

Valid_gender_samples <- Reduce(intersect, samples_list)


#--------

#-------- drop the gene and transcript identifiers and create a matrix

# the first column will be either transcript or column ids
feature_identifier <- colnames(ge)[1]

# drop any columns that contain other identifers
ge <- ge %>%
  dplyr::select(!!rlang::sym(feature_identifier), dplyr::starts_with("BS_")) %>%
  dplyr::filter(complete.cases(.)) %>%
  tibble::column_to_rownames(var = feature_identifier)

# create a matrix
gene_expression_mat <- as.matrix(ge)

# filter matrix columns to samples with valid values in all targetColumns
gene_expression_mat <- gene_expression_mat[ , colnames(gene_expression_mat) %in% Valid_gender_samples]
#--------

#--------Prepare input data matrix and response columns for model build

#need samples x transcripts for input to predictive model

input_mat <- t(gene_expression_mat)

#Extract target_column_list values from histologies.
#Put these values in a dataframe, targets_unsorted, with columns Kids_First_Biospecimen_ID, 
#target_column_list[[1]][1], target_column_list[[1]][2], ...
#So, we need a character vecor with the proposed column names for dataframe targets_unsorted

targets_unsorted_columns <- vector(mode="character", length=length(target_column_list[[1]]) + 1)
targets_unsorted_columns[1] <- "Kids_First_Biospecimen_ID"
for (t in seq(1:length(target_column_list[[1]]))) {
  targets_unsorted_columns[t + 1] <- target_column_list[[1]][t]
  
}

targets_unsorted <- histologies[histologies$Kids_First_Biospecimen_ID %in% rownames(input_mat),
                               targets_unsorted_columns]

#Match sequence of rownames(input_mat) and targets[, 1].
#targets holds the target_column_list values in the same Kids_First_Biospecimen_ID sequence as rownames(input_mat).

match_index <- unlist(sapply(rownames(input_mat), function(x) which(targets_unsorted[, 1] == x)))

targets <- targets_unsorted[match_index, c(1:ncol(targets_unsorted))]
#--------

#--------Partition input_mat and targets into train and test sets
# depending on train_percent argument generate four files: train_expression, test_expression, train_targets, test_targets
# or two files: train_expression, train_targets

if (train_percent < 1) {
  train_set_index <- sample(1:nrow(input_mat), floor(train_percent*nrow(input_mat)))
  test_set_index <- setdiff(1:nrow(input_mat), train_set_index)

  train_expression <- input_mat[train_set_index, ]
  test_expression <- input_mat[test_set_index, ]

  train_targets <- targets[train_set_index, ]
  test_targets <- targets[test_set_index, ]

  # Fill a column in targets that designates which set a sample is in (train vs. test).
  targets$train_or_test <- NA
  targets[train_set_index, "train_or_test"] <- "train"
  targets[test_set_index, "train_or_test"] <- "test"
} else {
  train_expression <- input_mat
  train_targets <- targets

  # all smaples are in train set here
  targets$train_or_test <- "train"
}
#--------

#--------Save train, (test, if used) and targets files to output directory

if (train_percent < 1) {
  train_file <- file.path(output_directory, opt$train_expression_file_name)
  test_file <- file.path(output_directory, opt$test_expression_file_name)
  saveRDS(train_expression, train_file)
  saveRDS(test_expression, test_file)
} else {
  train_file <- file.path(output_directory, opt$train_expression_file_name)
  saveRDS(train_expression, train_file)
}

if (train_percent < 1) {
  targets_file <- file.path(output_directory, opt$full_targets_file_name)
  train_targets_file <- file.path(output_directory, opt$train_targets_file_name)
  test_targets_file <- file.path(output_directory, opt$test_targets_file_name)

  write_tsv(targets, targets_file,
            na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(train_targets, train_targets_file,
            na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(test_targets, test_targets_file,
            na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
} else {
  targets_file <- file.path(output_directory, opt$full_targets_file_name)
  train_targets_file <- file.path(output_directory, opt$train_targets_file_name)

  write_tsv(targets, targets_file,
            na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(train_targets, train_targets_file,
            na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
}


cat(paste("\n\nEnd script 01_clean_split_data at", Sys.time(), "\n", sep=" "))
