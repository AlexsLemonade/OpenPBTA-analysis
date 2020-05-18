
#--------Library calls

library(glmnet)
library(glmnetUtils)
library(readr)
library(caret)

#--------

# This script is the third of a three-script workflow to predict sex from gene expression.
# We evaluate an elastic net logistic regression model that was trained in script 02-train_elasticnet.R
# from a training set and known target values set that was cleaned and saved in script 01-clean_split_data.R.

# Command-line arguments for this script identify the test set input directory where the test and target sets out of
# 01-clean_split_data.R are located, the model object input directory where the model object and transcript index files
# out of 02-train_elasticnet are located, the output directory where the evaluation files will be saved,
# the file name components (filename_lead, seed, MAD threshold) necessary to create
# the required file paths for input and output.

#________ Declare command line arguments

option_list <- list(
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory"
  ),
  optparse::make_option(
    c("-c", "--test_target_column"),
    type = "character",
    default = NULL,
    help = "column to compare for prediction accuracy"
  ),
  optparse::make_option(
    c("-n", "--test_expression_file_name"),
    type = "character",
    default = NULL,
    help = "input file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-r", "--test_targets_file_name"),
    type = "character",
    default = NULL,
    help = "input file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-m", "--model_object_file_name"),
    type = "character",
    default = NULL,
    help = "input file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-t", "--model_transcripts_file_name"),
    type = "character",
    default = NULL,
    help = "input file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-u", "--cm_set_file_name"),
    type = "character",
    default = NULL,
    help = "output file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-v", "--cm_file_name"),
    type = "character",
    default = NULL,
    help = "output file name constructed from arguments used in scripts 01 and 02"
  ),
  optparse::make_option(
    c("-w", "--summary_file_name"),
    type = "character",
    default = NULL,
    help = "output file name constructed from arguments used in scripts 01 and 02"
  )

)


# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

cat(paste("\n\nStart script 03-evaluate_model; output=", opt$output_directory, "at", Sys.time(), "\n", sep=" "))

# Create specified output directory if it does not yet exist
output_directory <- opt$output_directory
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

#--------

#--------test for the existence of the expected training and target value sets.

test_expression_file_name <- file.path(opt$test_expression_file_name)
test_targets_file_name <- file.path(opt$test_targets_file_name)

# print error and quit if not found.
if (!file.exists(test_expression_file_name)) {
  stop(paste("test set ", test_expression_file_name, "does not exist. Check all arguments."))
}

print("test set found!")

if (!file.exists(test_targets_file_name)) {
  stop(paste("targets set", test_targets_file_name, "does not exist. Check all arguments.", sep=" "))
}

print("targets set found!")

# safe to read if we get here

test_expression <- readRDS(test_expression_file_name)
test_targets <- read.delim(test_targets_file_name, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# check for presence of --test_target_column
if (!(opt$test_target_column %in% colnames(test_targets))) {
  stop(paste("target column ", opt$test_target_column, "does not exist. Check all arguments."))
}

print("target column found!")
#--------

#--------test for the existence of the expected model and transcript index files

model_file_name <- file.path(opt$model_object_file_name)
transcript_index_file_name <- file.path(opt$model_transcripts_file_name)

# print error and quit if not found.
if (!file.exists(model_file_name)) {
  stop(paste("model file ", model_file_name, "does not exist. Check all arguments."))
}

print("model file found!")

if (!file.exists(transcript_index_file_name)) {
  stop(paste("transcript index file", transcript_index_file_name, "does not exist. Check all arguments.", sep=" "))
}

print("transcript index set found!")

# safe to read if we get here

best_fit <- readRDS(model_file_name)
model_transcripts <- readRDS(transcript_index_file_name)

#--------

#-------- predict from model and evaluate sex predictions with caret::confusionMatrix

predict_class <- predict(best_fit, newx = test_expression[, model_transcripts], type = "class", s = best_fit$lambda.1se)
predict_prob <- predict(best_fit, newx = test_expression[, model_transcripts], type = "response", s = best_fit$lambda.1se)

cm_set <- data.frame(obs = test_targets[, opt$test_target_column])
cm_set$Male = unname(predict_prob)
cm_set$Female = 1 - unname(predict_prob)
cm_set$pred <- factor(unname(predict_class))
cm_set$sample <- test_targets[, "Kids_First_Biospecimen_ID"]

cm <- tryCatch(confusionMatrix(data = cm_set$pred, reference = cm_set$obs),
               error=function(e) e)
if (!inherits(cm, "error")) {
  cat("\n\n")
  print(cm)
  
  cm_set_file <- file.path(output_directory, opt$cm_set_file_name)
  
  write_tsv(cm_set, cm_set_file,
                na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  
  cm_file <- file.path(output_directory, opt$cm_file_name)
  saveRDS(cm, cm_file)
}


two_class_summary <- tryCatch(twoClassSummary(cm_set, lev = levels(cm_set$obs)),
                              error=function(e) e)
if (!inherits(two_class_summary, "error")) {
  cat("\n\n")
  print(two_class_summary)
  
  summary_file <- file.path(output_directory, opt$summary_file_name)
  saveRDS(two_class_summary, summary_file)
}

#--------

cat(paste("\n\nEnd script 03-evaluate_model; output=", opt$output_directory, "at", Sys.time(), "\n", sep=" "))
