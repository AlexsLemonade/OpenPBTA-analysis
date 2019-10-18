#Functionality/purpose:
  
#This will clean the gene expression data -- i.e., it drops anything with invalid labels 
#in either reported_gender or germline_sex_estimate and if the training partition size < 1.0, 
#it will split the data into training and testing sets that can be saved as separate files to be used downstream.

#Inputs and arguments:
#gene expression matrix
#histologies file
#seed with some default
#training partition size with some default
#output filenames or output directory, depending on how it is implemented


#--------Library calls and run parameters.  Modify the assignments as needed for subsequent runs.

library(dplyr)
library(tibble)
library(readr)

# We tune our predictive parameters using cross-validation on a training set and test the predictive accuracy
# of the optimal model on a holdout test set.
# train_percent represents the train/test partition; the seed is used to select training samples randomly.
# The test set is whatever did not get selected for the training set.

train_percent <- 0.70
set.seed(36354)

# Several gene expression files can serve as input to this analysis.
# Assign the file name for this run to gene_expression_file_name below.

gene_expression_file_name <- "pbta-gene-expression-kallisto.stranded.rds"

# training and test data files are saved into the directory named below.
# A TSV file that contains the biospecimen ID, reported_gender, and germline_sex_estimate 
# information that has been cleaned is also saved to this directory.

output_directory <- "scratch"
#--------

#--------Read input data files

#RNA-Seq data is in pbta-gene-expression-kallisto.rds.

ge <- readRDS(paste("../../data/", gene_expression_file_name, sep=""))


#Metadata is in pbta-histologies.tsv.  reported_gender is our target variable. 
#Kids_First_Biospecimen_ID is the unique identifier for each sample. 

histologies <- read.delim("../../data/pbta-histologies.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

valid_reported_gender_samples <- histologies[which(toupper(histologies[, "reported_gender"]) == "MALE" | 
              toupper(histologies[, "reported_gender"]) == "FEMALE"), "Kids_First_Biospecimen_ID"]

valid_germline_sex_estimate_samples <- histologies[which(toupper(histologies[, "germline_sex_estimate"]) == "MALE" | 
              toupper(histologies[, "germline_sex_estimate"]) == "FEMALE"), "Kids_First_Biospecimen_ID"]

Valid_gender_samples <- intersect(valid_reported_gender_samples, valid_germline_sex_estimate_samples)
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

# filter matrix columns to sample with valid reported_gender and germline_sex_estimate values
gene_expression_mat <- gene_expression_mat[ , colnames(gene_expression_mat) %in% Valid_gender_samples]
#--------

#--------Prepare input data matrix and response columns for model build

#need samples x transcripts for input to predictive model

input_mat <- t(gene_expression_mat)

#Extract reported_gender and germline_sex_estimate values from histologies.  
#Put these value in a three-column dataframe c("Kids_First_Biospecimen_ID", "reported_gender", "germline_sex_estimate").

targets_unsorted <- histologies[histologies$Kids_First_Biospecimen_ID %in% rownames(input_mat), 
                               c("Kids_First_Biospecimen_ID", "reported_gender", "germline_sex_estimate")]

#Check sequence of rownames(input_mat) and targets[, 1].  
#reported_gender_response holds the reported_gender values in the same Kids_First_Biospecimen_ID sequence as rownames(df).

match_index <- unlist(sapply(rownames(input_mat), function(x) which(targets_unsorted[, 1] == x)))

targets <- targets_unsorted[match_index, c(1:ncol(targets_unsorted))]
#--------

#--------Partition input_mat into train and test sets

if (train_percent < 1) {
  train_set_index <- sample(1:nrow(input_mat), floor(train_percent*nrow(input_mat)))
  test_set_index <- setdiff(1:nrow(input_mat), train_set_index)
  
  train_set <- input_mat[train_set_index, ]
  test_set <- input_mat[test_set_index, ]
} else {
  train_set <- input_mat
}
#--------

#--------Save train, (test, if used) and targets files to output directory

if (train_percent < 1) {
  save(train_set, test_set, file=paste("../../", output_directory, "/", "train_test_data.RData", sep=""))
} else {
  save(train_set, file=paste("../../", output_directory, "/", "train_data.RData", sep=""))
}

write_tsv(targets, paste("../../", output_directory, "/", "targets.tsv", sep=""),
          na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

