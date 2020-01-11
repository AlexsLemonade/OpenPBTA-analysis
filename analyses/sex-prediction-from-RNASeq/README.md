# Sex Prediction from RNASeq

The 01-clean_split_data.R, 02-train_elasticnet.R, 03-evaluate_model.R, 04-present_results.Rmd pipeline trains and evaluates an elasticnet logistic regression model to predict sex from RNASeq data.
The training features are gene expression transcripts, and the training labels are reported_gender values for each sample.

The pipeline is a response to issue [#84](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/84) which was raised to check whether, in some histologies, silencing might be breaking down, potentially resulting in changes in X inactivation.
Based on the accuracy achieved here, this is probably not happening, and this classifier can be helpful to predict values for datasets without annotated sex information.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  

- [How to run sex prediction from RNASeq](#how-to-run-sex-prediction-from-seq)
- [Summary of Methods](#summary-of-methods)
  - [Train/Test](#train-test)
  - [Model Building](#model-building)
  - [Model Evaluation](#model-evaluation)
- [General usage of scripts](#general-usage-of-scripts)
  - [run-sex-prediction-from-RNASeq.sh](#run-sex-prediction-from-seq)
  - [01-clean_split_data.R](#01-clean-split-data)
  - [02-train_elasticnet.R](#02-train_elasticnet)  
  - [03-evaluate_model.R](#03-evaluate_model)
  - [04-present_results.Rmd](#04-present_results)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## How to run sex-prediction-from-seq

To run the full pipeline of data prep, model building, model evaluation and results presentation, call the bash script:

```
bash run-sex-prediction-from-RNASeq.sh
```
This bash script requires arguments that are set in the USER-SPECIFIED ARGUMENTS section of the script.  
The script will return:

- Plots and tables in a notebook: `04-present_results.html`.  
  - Tables: confusion matrix and two class summary reports for the model with maximum accuracy.  
  - Plots: strength of calls on the test set for the maximum accuracy model, predictive accuracy vs. number of training transcripts and number non-zero features vs. number of training transcripts.
- Files organized in subfolders of the folder containing run-sex-prediction-from-RNASeq.sh as follows:
  - `processed_data` holds the output of script 01-clean_split_data.R.
  Five files: expression and target files for training and testing plus a file containing the training and testing target values combined.
  - `models` holds the output of script 02-train_elasticnet.R.
  Three files for each value in the run-sex-prediction-from-RNASeq.sh argument TRANSCRIPT_TAIL_PERCENT_ARRAY:
    - the best fitting model object for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value, 
    - the indices of the training transcripts used in that model
    - the non-zero coefficients of the model.
  - `results` holds the output of script 03-evaluate_model.R.  
  Three files for each value in the run-sex-prediction-from-RNASeq.sh argument TRANSCRIPT_TAIL_PERCENT_ARRAY: 
    - the caret::confusionMatrix for the model for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value, 
    - the caret::twoClassSummary object for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value (Note: caret::twoClassSummary function can fail under certain conditions, so twoClassSummary files may be missing for certain TRANSCRIPT_TAIL_PERCENT_ARRAY values.)
    - prediction probabilities for each sample in the test set.
  - `results/plots` holds .png copies of the plots presented in the `04-present_results.html` notebook.
  
## Summary of Methods

### train-test

Input data is split into training and test partitions according to the user-specified argument TRAIN_PERCENT.


### model-building

The glmnet package is used to fit an elastic net logistic regression model via penalized maximum likelihood.
The glmnetUtils package is used to do elastic net cross-validation for alpha and lambda simultaneously.

### model-evaluation

The caret package is used to generate a confusion matrix object and a two-class summary object for each model.
Definitions of the statistics presented are found in the documentation for the [confusion matrix function](http://topepo.github.io/caret/measuring-performance.html#measures-for-predicted-classes) and the [two-class summary function](http://topepo.github.io/caret/measuring-performance.html#measures-for-class-probabilities).


## General usage of scripts

**Overall notes about these scripts:**
- The scripts are sequential as noted by their number.
- All file path-related arguments assume the file path specified is relative to `OpenPBTA-analysis/analyses/sex-prediction-from-RNASeq`.
- The scripts create user-specified output directories that do not exist at run time.
- The scripts add files whose names do not match existing files, if any, in the output directories.  
Output files whose names match existing files overwrite the previous versions.

### run-sex-prediction-from-seq

A bash shell script that runs the entire pipeline.  
Global arguments for the pipeline are specified in the USER-SPECIFIED ARGUMENTS section of the script.


**Argument descriptions**
```
PROCESSED
output directory of script 01, input directory of scripts 02 and 03

MODELS
outtput directory of script 02, input directory of script 03

RESULTS
output directory of script 03

SEED
seed for R's pseudorandom number generator

FILENAME_LEAD
argument for script 01 output file specification and script 02 and 03 input and output file specification

TRANSCRIPT_TAIL_PERCENT_ARRAY
a bash shell script array to control median absolute deviation filtering of training transcripts.
Array values must be between 0 and 1.
No limit on the size of the array.
For each value in the array, a model is trained and evaluated by filtering out the bottom (1 - TRANSCRIPT_TAIL_PERCENT_ARRAY value) of transcripts from the training set.

TRAIN_PERCENT
a value greater than 0 and less than or equal to 1.
The expression dataset is partitioned TRAIN_PERCENT for training and (1 - TRAIN_PERCENT) for testing.
If TRAIN_PERCENT = 1, no test set is created.
  
TRAIN_TARGET_COLUMN
argument to specify what column in the target data frame will be used as labels during training.

targetColumns
a bash shell script array specifying the columns to be used as ground truth labels when evaluating a model's 
performance.
Performance reports are produced for each value in the array.  No limit on the size of the array.

```

### 01-clean-split-data

This script cleans the gene expression data -- i.e., it drops anything with invalid labels in either 
reported_gender or germline_sex_estimate, and if the training partition size < 1.0, it splits the data 
into training and testing sets that are saved as separate files to be used downstream.


**Argument descriptions**
```
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds 
    File path and file name of the OpenPBTA gene expression file used for the current run of the pipeline(RDS)
  --metadata ../../data/pbta-histologies.tsv 
    File path and file name of the OpenPBTA histologies file used for the current run(TSV)
    
Values for the following arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS.

  --output_directory PROCESSED 
    File path where you would like the cleaned training and test sets to be stored.
  --train_expression_file_name TRAIN_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_expression.RDS 
    Name for the training expression set output file
  --test_expression_file_name TEST_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_expression.RDS
    Name for the test expression set output file
  --train_targets_file_name TRAIN_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_targets.tsv
    Name for the training label set output file
  --test_targets_file_name TEST_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_targets.tsv
    Name for the test label set output file
  --full_targets_file_name FULL_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_full_targets.tsv
    Name for the combined training and test label set output file
  --seed $SEED
    seed for random number generation
  --train_percent $TRAIN_PERCENT 
    controls training/test partition

```

### 02-train_elasticnet

This script builds elastic net logistic regression models from the training set and known target values set
that were cleaned and saved in script 01-clean_split_data.R.
The models are generated by a loop over the values of the TRANSCRIPT_TAIL_ARRAY argument in the USER_ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R.  
i is the index of the loop over the values of TRANSCRIPT_TAIL_ARRAY.

 --train_expression_file_name ${PROCESSED}/$TRAIN_EXPRESSION_FILE_NAME 
 --train_targets_file_name ${PROCESSED}/$TRAIN_TARGETS_FILE_NAME 
 --output_directory $MODELS 
 --model_object_file_name $MODEL_OBJECT_FILE_NAME 
 --model_transcripts_file_name $MODEL_TRANSCRIPTS_FILE_NAME 
 --model_coefs_file_name $MODEL_COEFS_FILE_NAME 
 --train_target_column $TRAIN_TARGET_COLUMN 
 --transcript_tail_percent ${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}
 
```
### 03-evaluate_model

This script evaluates the elastic net logistic regression models that were trained in script 
02-train_elasticnet.R from a training set and known target values set that was cleaned and saved 
in script 01-clean_split_data.R.
This script runs within the same loop described for script 02-train_elasticnet.R, and within a second nested loop over the values of the targetColumns argument in the USER_SPECIFIED_ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.
Output files for each value of targetColumns go to a separate directory named with the current targetColumns value.

Model evaluations are done with the caret::confusionMatrix and caret::twoClassSummary functions.
twoClassSummary can fail if one of the Confusion Matrix rows is all zeros.
In that case, no twoClassSummary file is output for the transcript_tail_percent value in question.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R and 02-train_elasticnet.R.
i is the index of the loop over the values of TRANSCRIPT_TAIL_ARRAY.
t is the index of the loop over the values of targetColumns.

      --test_expression_file_name ${PROCESSED}/$TEST_EXPRESSION_FILE_NAME 
      --test_targets_file_name ${PROCESSED}/$TEST_TARGETS_FILE_NAME 
      --model_object_file_name ${MODELS}/$MODEL_OBJECT_FILE_NAME 
      --model_transcripts_file_name ${MODELS}/$MODEL_TRANSCRIPTS_FILE_NAME 
      --test_target_column ${t} 
      --output_directory $RESULTS_OUTPUT_DIRECTORY=${RESULTS}/${TEST_TARGET_COLUMN}
      --cm_set_file_name $CM_SET_FILE=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_prediction_details.tsv
        A file containing prediction proabilities for each sample in the test set.
      --cm_file_name $CM_SET=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_confusion_matrix.RDS
        A file containing the confusion matrix object out of caret::confusionMatrix.
      --summary_file_name $SUMMARY_FILE=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_two_class_summary.RDS
        A file containing the two class summary object out of caret::twoClassSummary.

```

### 04-present_results

A Notebook that produces tables and plots from the output files of 02-train_elasticnet.R and 03-evaluate_model.R.  
Plots produced are Strength of Calls at Maximum Accuracy, Predictive Accuracy vs. Number of Training Transcripts,
Number of Non-Zero Features vs. Number of Training Transcripts.
Tables are Confusion Matrix at Maximum Accuracy and Two Class Summary at Maximum Accuracy.
Each of the tables and charts is produced for each value of the targetColumns array specified in the USER-SPECIFIED ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R, 02-train_elasticnet.R
and 03-evaluate_model.R.

  --results_dir=${RESULTS} 
  --model_dir=$MODELS 
  --cm_set=$CM_SET 
  --seed=$SEED 
  --target_columns=$targetColumns_to_pass
    A string constructed from the values of the targetColumns array in USER-SPECIFIED ARGUMENTS.

```
