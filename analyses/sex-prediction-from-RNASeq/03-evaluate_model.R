#--------Library calls

cat(paste("\n\nStart script 03-evaluate_model at", Sys.time(), "\n", sep=" "))

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
    c("-i", "--test_set_input_directory"),
    type = "character",
    default = NULL,
    help = "input directory"
  ),
  optparse::make_option(
    c("-j", "--model_input_directory"),
    type = "character",
    default = NULL,
    help = "input directory"
  ),  
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory"
  ),
  optparse::make_option(
    c("-f", "--filename_lead"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  ),
  optparse::make_option(
    c("-c", "--target_column"),
    type = "character",
    default = "reported_gender",
    help = "A character vector identifying column being predicted"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 36354,
    help = "seed integer",
    metavar = "integer"
  ),  
  optparse::make_option(
    c("-p", "--transcript_tail_percent"),
    type = "double",
    default = 0.25,
    help = "Filter out the bottom (1 - transcript_tail_percent) of transcripts from the training set",
    metavar = "double"
  )
  
)


# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Create specified output directory if it does not yet exist
output_directory <- opt$output_directory
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

#--------

#--------test for the existence of the expected training and target value sets. 

test_expression_file_name <- file.path(opt$test_set_input_directory, paste(opt$filename_lead, opt$seed, "test_expression.RDS", sep = "_"))
test_targets_file_name <- file.path(opt$test_set_input_directory, paste(opt$filename_lead, opt$seed, "test_targets.tsv", sep = "_"))

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

# check for presence of --target_column
if (!(opt$target_column %in% colnames(test_targets))) {
  stop(paste("target column ", opt$target_column, "does not exist. Check all arguments."))
}

print("target column found!")
#--------

#--------test for the existence of the expected model and transcript index files

model_file_name <- file.path(opt$model_input_directory, paste(opt$filename_lead, opt$seed, opt$transcript_tail_percent, "model_object.RDS", sep = "_"))
transcript_index_file_name <- file.path(opt$model_input_directory, paste(opt$filename_lead, opt$seed, opt$transcript_tail_percent, "model_transcripts.RDS", sep = "_"))

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

cm_set <- data.frame(obs = test_targets[, opt$target_column])
cm_set$Male = unname(predict_prob)
cm_set$Female = 1 - unname(predict_prob)
cm_set$pred <- factor(unname(predict_class))
cm_set$sample <- test_targets[, "Kids_First_Biospecimen_ID"]

cm <- confusionMatrix(data = cm_set$pred, reference = cm_set$obs)
cat("\n\n")
cm

two_class_summary <- twoClassSummary(cm_set, lev = levels(cm_set$obs))
cat("\n\n")
two_class_summary

cm_set_file <- file.path(opt$output_directory, paste(opt$filename_lead, opt$seed, opt$transcript_tail_percent, 
                                                 "prediction_details.tsv", sep = "_"))
write_tsv(cm_set, cm_set_file,
          na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")


cm_file <- file.path(opt$output_directory, paste(opt$filename_lead, opt$seed, opt$transcript_tail_percent, 
                                                          "confusion_matrix.RDS", sep = "_"))
saveRDS(cm, cm_file)

summary_file <- file.path(opt$output_directory, paste(opt$filename_lead, opt$seed, opt$transcript_tail_percent, 
                                                      "two_class_summary.RDS", sep = "_"))
saveRDS(two_class_summary, summary_file)


#--------

cat(paste("\n\nEnd script 03-evaluate_model at", Sys.time(), "\n", sep=" "))